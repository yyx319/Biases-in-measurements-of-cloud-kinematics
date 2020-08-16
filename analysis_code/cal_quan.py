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
mH = 2.34e-24
from scipy.interpolate import RegularGridInterpolator
from yt.units.yt_array import YTQuantity as q
from scipy import stats
from scipy import interpolate
pc = 3.0856e18
km = 1.e5
cm_to_pc = q(1.0/pc, 'pc/cm')
cm_to_km = q(1.0/km, 'km/cm')
kmspc_to_pers = q(km/pc, 'pc/km')
dvdr_unit = q(1.0, 'km/s/pc')
# scale parameter den = den*sca       
tau=0
dens=0
ACF=1
# simulation parameter
beta_a = ['0.2','2','20']
time_a = ['10','30','60']

directory = '/home/yuxuan/tracer_ls_rel'
#setting LOS
prj_ori_a = ['x','y','z'] #'x','y','z'
sp1_a = ['CO','CO_4_3','C18O','pNH3','HCN','C18O_4_3','N2H']
sp2_a = ['co','co_4_3','c18o','pnh3','hcn','c18o_4_3','n2h']

tau_m = np.zeros((3,3,3,7))
rho_L = np.zeros((3,3,3,7))


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
    
for i, beta in enumerate(beta_a):
    for j,time in enumerate(time_a):
        exec(open("/avatar/yuxuan/tracer_ls_rel/C12_DATA/ReadData.py").read())
        for k,prj_ori in enumerate(prj_ori_a):  
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
            
            for s, sp1, sp2 in zip(range(len(sp1_a)), sp1_a, sp2_a):
                if tau == 1:
                    exec(open("/avatar/yuxuan/tracer_ls_rel/molecule/tau_field.py").read())
                    field = sp1+"_tau" 
                    weight =  "density" 
                    
                    ad = ds.all_data()  # This is a region describing the entire box,
                                    # but note it doesn't read anything in yet!
                    t_m = ad.quantities.weighted_average_quantity(field, weight)
                    tau_m[i,j,k,s]=t_m.value  
                    
                if dens ==1 or ACF==1:
                    exec(open("/avatar/yuxuan/tracer_ls_rel/molecule/lum_field.py").read())
                    ad = ds.all_data()
                    
                if dens ==1:
                    field = "density" 
                    weight = sp1+"_Luminosity_"+prj_ori 
                    r_L = ad.quantities.weighted_average_quantity(field, weight)
                    r_L = r_L.value/mH
                    rho_L[i,j,k,s]=r_L  
                    
                if ACF == 1:
                    locals()['L_%d_%d_%d_%d'%(i,j,k,s)] = ad[sp1+'_Luminosity_'+prj_ori]
                    locals()['L_%d_%d_%d_%d'%(i,j,k,s)] = np.reshape(locals()['L_%d_%d_%d_%d'%(i,j,k,s)],(256,256,256))

if tau==1:
    print(tau_m)
    np.save("%s/tau_m"%(directory), tau_m)
if dens==1:
    print(rho_L)
    np.save("%s/rho_L"%(directory), rho_L) 

if ACF==1:
    ACL = np.zeros((3,3,5))
    for i, beta in enumerate(beta_a):
        for j,time in enumerate(time_a): 
            for s, sp1, sp2 in zip(range(len(sp1_a)), sp1_a, sp2_a):
                L = ( locals()['L_%d_%d_0_%d'%(i,j,s)]+locals()['L_%d_%d_1_%d'%(i,j,s)]+locals()['L_%d_%d_2_%d'%(i,j,s)] )/3
                fL = np.fft.fftn(L, s=None, axes=None)
                psi = fL*np.conj(fL)
                A3 = np.fft.ifftn(psi)
                A3 = np.real(A3)/256**3
                A3_nor = A3/A3[0,0,0]
                A3_nor = np.roll(A3_nor, [128,128,128], axis=[0,1,2])
                dx = 4.6/256
                x = dx * (np.arange(256) - 128)          # Get x position in rolled A3 array; here dx is the size of a cell
                xxx, yyy, zzz = np.meshgrid(x,x,x)      # Get x, y, z positions of every cell in the A3 array
                r = np.sqrt(xxx**2 + yyy**2 + zzz**2)  # Get radial distance of every point in the A3 array
                A1_nor, r_bins, bin_numbers = stats.binned_statistic(np.ravel(r), np.ravel(A3_nor), bins=224)    # Compute mean value of A3_nor in bins of radius
                lag = 0.5 * (r_bins[1:] + r_bins[:-1])
                np.savetxt('%s/ACF/beta%s_%s/A1_%s.txt'%(directory, beta, time, sp2), A1_nor)
                    
    for i, beta in enumerate(beta_a):
        for j,time in enumerate(time_a): 
            for s, sp1, sp2 in zip(range(len(sp1_a)), sp1_a, sp2_a):
                A1_nor = np.loadtxt('%s/ACF/beta%s_%s/A1_%s.txt'%(directory, beta, time, sp2))
                f = interpolate.interp1d(A1_nor, lag)
                if s!=0 and s!=1:
                    acl = f(0.5)
                    ACL[i,j,s-1] = acl
            
    print(ACL)
    np.save("%s/ACL"%(directory), ACL) 


