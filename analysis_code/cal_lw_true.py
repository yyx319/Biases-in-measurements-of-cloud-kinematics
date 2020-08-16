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
beta_a = ['0.2','2','20']
time_a = ['10','30','60']

#setting LOS
prj_ori_a = ['x','y','z'] #'x','y','z'
theta_a = np.array([90., 90., 0.])*np.pi/180. #90., 90., 0.
phi_a = np.array([0., 90., 0.])*np.pi/180. #0., 90., 0.

prj_ori_a = prj_ori_a[2:3]
theta_a = theta_a[2:3]
phi_a = phi_a[2:3] 



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
             
        for prj_ori, theta, phi in zip(prj_ori_a, theta_a, phi_a):
            LOS = [np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(phi)]
            ####################
            #HCN_Luminosity PPV#
            ####################
            from yt.config import ytcfg
            from yt.analysis_modules.ppv_cube.api import PPVCube
            import yt.units as u
            cube = PPVCube(ds, LOS, 'density', (-vr,vr,res_v,"km/s"), dims=256, method="sum")
            cube.write_fits("/avatar/yuxuan/tracer_ls_rel/PPV/beta%s_%s_%s/cube_%s_%s_%s_dens.fits"%(beta,time,prj_ori,beta, time, prj_ori), overwrite=True, length_unit="pc")
            
            from astropy.io import fits
            cube = fits.open("/avatar/yuxuan/tracer_ls_rel/PPV/beta%s_%s_%s/cube_%s_%s_%s_dens.fits"%(beta,time,prj_ori,beta, time, prj_ori))   
            cube = cube['density'].data
            cube = cube[:,:-1,:-1]
            directory ='/avatar/yuxuan/tracer_ls_rel/moment'
            v = np.linspace(-4, 4, 201)
            v = (v[1:]+v[:-1])/2
            v = v[:,None, None]
            m0 = np.sum(cube, axis=0)
            np.savetxt("%s/beta%s_%s_%s/m0_true.txt"%(directory,beta,time,prj_ori), m0)
            
            m1 = np.sum(cube*v, axis=0)/m0
            np.savetxt("%s/beta%s_%s_%s/m1_true.txt"%(directory,beta,time,prj_ori), m1)
           
            m2 = np.sum((v-m1)**2*cube, axis=0)/m0
            m2 = np.sqrt(m2)  
            np.savetxt("%s/beta%s_%s_%s/m2_true.txt"%(directory,beta,time,prj_ori), m2)
                 

