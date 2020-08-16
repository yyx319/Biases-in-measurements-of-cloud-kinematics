# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 09:40:39 2019

@author: YYX
"""
########
#set up#
########
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
from scipy.interpolate import RegularGridInterpolator
plt.switch_backend('agg')
import h5py
import yt
from yt import derived_field
from scipy.interpolate import RectBivariateSpline as rbs
from yt.units.yt_array import YTQuantity as q

pi = np.pi
km = 1.e5
pc = 3.0856e18
yr = 365.25*24.*3600.

c = 2.99792458e10
kB = 1.38065812e-16
mH = 2.34e-24	#effective mass per H nucleus
nH = 1.0e4  	#threshold volume density
t_rho = nH*mH		#threshold volume density in g/cm**3

unitconversion1 = q(1.0, 'K*km*pc*pc/erg')
unitconversion2 = q(1.0, 'erg/s/g')

cm_to_pc = q(1.0/pc, 'pc/cm')
cm_to_km = q(1.0/km, 'km/cm')
kmspc_to_pers = q(km/pc, 'pc/km')
dvdr_unit = q(1.0, 'km/s/pc')

# PPV parameter
vr=4 #velocity range -vr~+vr
res_v=200

# simulation parameter
beta_a = ['20']
time_a = ['60']
sca = 10
#setting LOS
prj_ori_a = ['x','y','z'] #'x','y','z'
theta_a = np.array([90., 90., 0.])*np.pi/180. #90., 90., 0.
phi_a = np.array([0., 90., 0.])*np.pi/180. #0., 90., 0.

prj_ori_a = prj_ori_a[0:1] #'x','z'
theta_a = theta_a[0:1]
phi_a = phi_a[0:1] 

sp1 = ['co','c18o','pnh3','hcn','n2h','co_4_3','c18o_4_3']
sp2 = ['CO','C18O','pNH3','HCN','N2H','CO_4_3','C18O_4_3']
abun = ['std','std','std','std','std','std','std']

sp_index_a = range(2,7)

second_m = []



# rebin the original PPV
def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = np.array(np.asarray(shape)/np.asarray(args), dtype=int)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)] 
    print(''.join(evList))
    return eval(''.join(evList))






for beta in beta_a:
    for time in time_a:
        ###########
        #READ Data#
        ###########
        # MRK
        # Put number of output file you want to read here. For example, to
        # read the file C12_Beta0.2_256_0000.h5, you put 00 here. To read the
        # file C12_Beta0.2_256_0010.h5, you put 10 here.
        frames= [int(time)] #[10,30,60]
        
        field_units = {'Bx':'code_magnetic','By':'code_magnetic','Bz':'code_magnetic',\
                       'velocity_x':'code_velocity','velocity_y':'code_velocity','velocity_z':'code_velocity',\
                       'density':'code_density'}
        unit_conversions = {} #to be read from the files disk
        # MRK
        # There are two things for you to fill in here.
        # directory = directory where the file you are trying to open is
        # located
        # setname_base = base name of file; for example, if you are trying
        # to open C12_Beta0.2_256_0000.h5, your setname_base =
        # 'C12_Beta0.2_256'
        directory = '/avatar/yuxuan/tracer_ls_rel/C12_DATA'
        setname_base = 'C12_Beta%s_256'%(beta)
        data_shape =np.array([256,256,256])
        
        for frame in frames:
            set_name = "%s/%s_%04d.h5"%(directory,setname_base,frame)
            fptr = h5py.File(set_name,"r")
            data={}
            try:
        
                print("reading data for %s"%set_name)
                for field in field_units:
                    data[field]=(fptr[field][:] ,field_units[field])
                    these_units = fptr[field+"_units"].value.decode('ascii').split(" ")
                    unit_conversions[field]=(float(these_units[0]),these_units[1])
                these_units = fptr['length_units'].value.decode('ascii').split(" ")
                unit_conversions['length']=(float(these_units[0]),these_units[1])
                bbox = np.array([[0.,1.],[0.,1.],[0.,1.]])
                ds = yt.load_uniform_grid(data,data_shape, bbox=bbox,periodicity=(True,True,True),\
                                          length_unit=unit_conversions['length'],\
                                          velocity_unit=unit_conversions['velocity_x'],\
                                          magnetic_unit=unit_conversions['Bx'],\
                                          mass_unit=unit_conversions['density'][0]*unit_conversions['length'][0]**3)
            except:
                raise
            finally:
                fptr.close()
            print("Computing plasma beta:")
            ad=ds.all_data()
            ad.get_data()  
            bx_mean_yt = np.mean(ad['Bx']).in_units('code_magnetic')
            by_mean_yt = np.mean(ad['By']).in_units('code_magnetic')
            bz_mean_yt = np.mean(ad['Bz']).in_units('code_magnetic')
            density = np.mean(ad['density'])
            cs = yt.YTArray(unit_conversions['velocity_x'][0],'cm/s')
            beta_ex = 8*np.pi*cs**2*density/(bx_mean_yt.in_units('gauss'))**2
            print(beta_ex)
            print(ad['density'])

            
        
        #### Units and constants #### 
       
        for prj_ori, theta, phi in zip(prj_ori_a, theta_a, phi_a):
            # add velocity gradient by smoothing over 
            res = 256
            sm_fac = 8 # smooth factor
            sm_res = int(res/sm_fac) # resolution of the smoothed cube         
            v = np.array(ad['velocity_'+prj_ori].reshape(res,res,res)*cm_to_km) # in km/s 
            x = np.array(ad['x']*cm_to_pc) # in pc
            y = np.array(ad['y']*cm_to_pc) 
            z = np.array(ad['z']*cm_to_pc) 
            xx = np.array(x.reshape(res,res,res)[:,0,0]) # in pc
            xx_sm = rebin(xx, sm_res)
            yy = np.array(y.reshape(res,res,res)[0,:,0]) # in pc
            yy_sm = rebin(xx, sm_res)
            zz = np.array(z.reshape(res,res,res)[0,0,:]) # in pc
            zz_sm = rebin(xx, sm_res)         
            dl = ad['dx'][0]*cm_to_pc # in pc
            dl_sm = sm_fac*dl 
            if prj_ori=='x':
                v_sm = rebin(v, sm_res, res, res) 
                v_grad_sm = np.zeros(v_sm.shape)
                v_grad_sm[1:-1,:,:] = (v_sm[2:,:,:] - v_sm[:-2,:,:]) / (2*dl_sm) # in km/s/pc
                v_grad_sm[0,:,:] = (v_sm[1,:,:] - v_sm[-1,:,:]) / (2*dl_sm)
                v_grad_sm[-1,:,:] = (v_sm[0,:,:] - v_sm[-2,:,:]) / (2*dl_sm)
                vgrad_interp_func = RegularGridInterpolator((xx_sm, yy, zz), v_grad_sm, bounds_error=False, fill_value=None)
            if prj_ori=='y':
                v_sm = rebin(v, res, sm_res, res) 
                v_grad_sm = np.zeros(v_sm.shape)
                v_grad_sm[:,1:-1,:] = (v_sm[:,2:,:] - v_sm[:,:-2,:]) / (2*dl_sm) # in km/s/pc
                v_grad_sm[:,0,:] = (v_sm[:,1,:] - v_sm[:,-1,:]) / (2*dl_sm)
                v_grad_sm[:,-1,:] = (v_sm[:,0,:] - v_sm[:,-2,:]) / (2*dl_sm)
                vgrad_interp_func = RegularGridInterpolator((xx, yy_sm, zz), v_grad_sm, bounds_error=False, fill_value=None)
            if prj_ori=='z':
                v_sm = rebin(v, res, res, sm_res) 
                v_grad_sm = np.zeros(v_sm.shape)
                v_grad_sm[:,:,1:-1] = (v_sm[:,:,2:] - v_sm[:,:,:-2]) / (2*dl_sm) # in km/s/pc
                v_grad_sm[:,:,0] = (v_sm[:,:,1] - v_sm[:,:,-1]) / (2*dl_sm)
                v_grad_sm[:,:,-1] = (v_sm[:,:,0] - v_sm[:,:,-2]) / (2*dl_sm)
                vgrad_interp_func = RegularGridInterpolator((xx, yy, zz_sm), v_grad_sm, bounds_error=False, fill_value=None)
            


            def _vgrad(field, data):
                x = np.ravel(data['x']*cm_to_pc) # in pc
                y = np.ravel(data['y']*cm_to_pc)
                z = np.ravel(data['z']*cm_to_pc) 
                grid = np.vstack((x, y, z)).T # array of the (x,y,z)
                vgrad = vgrad_interp_func( grid )
                vgrad = np.abs(vgrad)
                return vgrad.reshape(x.shape)*q(km/pc, '1/s')
                
            ds.add_field(('gas','vgrad_%s'%prj_ori), function=_vgrad, units="1/s")
            ds.periodicity = (True, True, True)
            
            
            for s in sp_index_a:
                lw_field = "%s_Luminosity"%(sp2[s])
                            
                #### Defining the field ####    
                file1 = "/avatar/yuxuan/tracer_ls_rel/molecule/n.txt"										#list of number densities in H/cm**3 (log scale)
                number = np.loadtxt(file1, delimiter=" ", dtype=None)
                file2 = "/avatar/yuxuan/tracer_ls_rel/molecule/dvdr.txt"										#list of velocity gradients in km/s/pc
                dvdr = np.loadtxt(file2, delimiter=" ", dtype=None)
                file3 = "/avatar/yuxuan/tracer_ls_rel/molecule/%s/%s_lum_std.txt"%(sp2[s], sp1[s])										#list of corresponding HCN luminosities per H molecule (log scale)
                lum_tab = np.loadtxt(file3, delimiter=" ", dtype=None)
                
                x = number[:]+np.log10(mH) 		#log scale of mass density in g/cm^3
                y = dvdr[:]-np.log10(pc/km) 	#log scale of velocity gradient in s^-1
                z = lum_tab[:,:] 					#log scale of luminosity produced by HCN per molecule of H in erg/s/H
                f = rbs(x,y,z)					#cubic spline interpolation
                
                def _Lum(field, data):
                    dvdr3d = data['vgrad_%s'%prj_ori]*np.sqrt(sca)
                    dens3d = data['density']*sca
                    mass = data['cell_mass']/np.sqrt(sca)
                    dvdr1d = np.ravel(dvdr3d)
                    dens1d = np.ravel(dens3d)
                    logdvdr = np.log10(dvdr1d)
                    logdens = np.log10(dens1d)
                    loglum = f(logdens,logdvdr, grid = False)
                    lum1d = 10**(loglum)
                    lum3d = lum1d.reshape(dens3d.shape)
                    l = mass*lum3d/mH
                    return l*unitconversion2
                
                	
                ds.add_field(('gas',lw_field), function=_Lum, units="erg/s", force_override=True)
                ds.periodicity = (True, True, True)
            
                LOS = [np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(phi)]
                ####################
                #HCN_Luminosity PPV#
                ####################
                from yt.config import ytcfg
                from yt.analysis_modules.ppv_cube.api import PPVCube
                import yt.units as u
                cube = PPVCube(ds, LOS, lw_field, (-vr,vr,res_v,"km/s"), dims=256, method="sum")
                cube.write_fits("/avatar/yuxuan/tracer_ls_rel/rescale_PPV/beta%s_%s_%s/cube_%.1f_%s_%s_%s_L%s.fits"%(beta,time,prj_ori,sca,beta, time, prj_ori, sp1[s]), overwrite=True, length_unit="pc")
                
                from astropy.io import fits
                cube = fits.open("/avatar/yuxuan/tracer_ls_rel/rescale_PPV/beta%s_%s_%s/cube_%.1f_%s_%s_%s_L%s.fits"%(beta,time,prj_ori,sca,beta,time,prj_ori,sp1[s]))   
                cube = cube['%s_Luminosity'%sp1[s]].data
                cube = cube[:,:-1,:-1]
                directory ='/avatar/yuxuan/tracer_ls_rel/rescale_moment'
                v = np.linspace(-4, 4, 201)
                v = (v[1:]+v[:-1])/2
                v = v[:,None, None]
                m0 = np.sum(cube, axis=0)
                np.savetxt("%s/beta%s_%s_%s/m0_%.1f_%s.txt"%(directory,beta,time,prj_ori,sca,sp1[s]), m0)
                
                m1 = np.sum(cube*v, axis=0)/m0
                np.savetxt("%s/beta%s_%s_%s/m1_%.1f_%s.txt"%(directory,beta,time,prj_ori,sca,sp1[s]), m1)
               
                m2 = np.sum((v-m1)**2*cube, axis=0)/m0
                m2 = np.sqrt(m2)  
                np.savetxt("%s/beta%s_%s_%s/m2_%.1f_%s.txt"%(directory,beta,time,prj_ori,sca,sp1[s]), m2)
                
                print(r'\beta=%s, t=%s t_{ff} orientation: %s'%(beta, time, prj_ori))
                lw = np.average(m2, weights=m0)
                second_m.append(lw)  
            np.savetxt("%s/linewidth_%.1f_%s_%s_%s.txt"%(directory,sca,beta,time, prj_ori), second_m)   


