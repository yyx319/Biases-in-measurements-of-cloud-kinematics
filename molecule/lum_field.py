# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 23:27:39 2019

@author: YYX
"""
import matplotlib.pyplot as plt
from matplotlib import pylab
import yt
from yt import derived_field
from scipy.interpolate import RectBivariateSpline as rbs
from scipy.interpolate import RegularGridInterpolator
from yt.units.yt_array import YTQuantity as q
import numpy as np
import h5py

#### Units and constants ####
# read data ds ad, set species sp1 sp2 and project orientation
pi = np.pi
km = 1.e5
pc = 3.0856e18
yr = 365.25*24.*3600.
c = 2.99792458e10
kB = 1.38065812e-16
mH = 2.34e-24	#effective mass per H nucleus
nH = 1.0e4  	#threshold volume density
t_rho = nH*mH		#threshold volume density in g/cm**3
nu_dict = {'CO':115.3e9,'C18O':109.78e9,'pNH3':23.69e9,'HCN':88.631847e9,'N2H':93.173809e9,'CO_4_3':461.0407682e9,'C18O_4_3':439.0887658e9}
nu = nu_dict[sp1]
lambd = c/nu	#wavelenght of HCN(1-0)
unitconversion1 = q(1.0, 'K*km*pc*pc/erg')
unitconversion2 = q(1.0, 'erg/s/g')

cm_to_pc = q(1.0/pc, 'pc/cm')
cm_to_km = q(1.0/km, 'km/cm')
kmspc_to_pers = q(km/pc, 'pc/km')
dvdr_unit = q(1.0, 'km/s/pc')

conversion = unitconversion1*lambd*lambd*lambd/8/pi/kB/pc/pc/km 	#convert erg/s to K km/s pc^2
locals()['conversion_'+sp1] = lambd*lambd*lambd/8/pi/kB/pc/pc/km 
bins = 100

#### Defining the field ####

file1 = "C:/Users/YYX/Desktop/Tracers_LS_relation/DESPOTIC/molecule/n.txt"								#list of number densities in H/cm**3 (log scale)
number = np.loadtxt(file1, delimiter=" ", dtype=None)
file2 = "C:/Users/YYX/Desktop/Tracers_LS_relation/DESPOTIC/molecule/dvdr.txt"									#list of velocity gradients in km/s/pc
dvdr = np.loadtxt(file2, delimiter=" ", dtype=None)
file3 = "C:/Users/YYX/Desktop/Tracers_LS_relation/DESPOTIC/molecule/%s/%s_lum_std.txt"%(sp1, sp2)									#list of corresponding HCN luminosities per H molecule (log scale)
Lum_tab = np.loadtxt(file3, delimiter=" ", dtype=None)

x = number[:]+np.log10(mH) 		# log scale of mass density in g/cm^3
y = dvdr[:]-np.log10(pc/km) 	# log scale of velocity gradient in s^-1
z = Lum_tab[:,:] 					# log scale of luminosity produced by molecule per molecule of H in erg/s/H
f = rbs(x,y,z)					# cubic spline interpolation


def _Lum(field, data):
    dvdr3d = data['vgrad_'+prj_ori]
    dens3d = data['density']
    mass = data['cell_mass']
    dvdr1d = np.ravel(dvdr3d)
    dens1d = np.ravel(dens3d)
    logdvdr = np.log10(dvdr1d)
    logdens = np.log10(dens1d)
    loglum = f(logdens,logdvdr, grid = False)
    lum1d = 10**(loglum)
    lum3d = lum1d.reshape(dens3d.shape)
    lum_sp = mass*lum3d/mH
    return lum_sp*unitconversion2
	
ds.add_field(('gas','%s_Luminosity'%sp1), function=_Lum, units="erg/s")
ds.periodicity = (True, True, True)




#### Density Profiles ####

rhomin,rhomax = ad.quantities.extrema("density")/(mH*q(1.0, 'g'))
mtotal = ad.quantities.total_quantity("cell_mass")
Lumtotal = ad.quantities.total_quantity("%s_Luminosity"%sp1)	#already defined
logmin = np.log10(rhomin)
logmax = np.log10(rhomax)
rhoarray = np.logspace(logmin,logmax,bins)
dlog = np.log10(rhomax/rhomin)/(bins-1)

massbin = []
locals()['%sbin'%sp1] = []
for i in range(0,bins,1):
    rho = ad["density"]/(mH*q(1.0, 'g'))
    if i==0:
        indices = np.where(rho <= rhoarray[i])
    if i==bins:
        indices = np.where(rho >= rhoarray[i])
    else:
        indices = np.where(np.logical_and(rho <= rhoarray[i], rho >= rhoarray[i-1])) 
    masssum = np.sum(ad["cell_mass"][indices])/(mtotal*dlog)
    lumsum = np.sum(ad['%s_Luminosity'%sp1][indices])/(Lumtotal*dlog)
    massbin.append(masssum)
    locals()['%sbin'%sp1].append(lumsum)