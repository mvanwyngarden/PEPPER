# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
from bisect import bisect_left
import pdb

def weights(Xiron):
    
    Xirons = np.arange(0, 1.025, 0.025)
    #returns the position within the list at which the given Xiron value would be inserted
    #if Xiron is equal to value in list it returns the index of the matching value
    pos = bisect_left(Xirons, Xiron)
    diff1 = abs(Xirons[pos]-Xiron)
    diff2 = abs(Xirons[pos-1]-Xiron)
    
    #if the given Xiron is equal to a value in the list the coefficient is 1
    if diff1 == 0: 
        coeff1 = 0
        coeff2 = 1
    else: 
        coeff1 = diff1/.025
        coeff2 = diff2/.025
  
    return pos, coeff1, coeff2
    
    
def tables():
    radtable = np.loadtxt(r'C:/Users/madir/Radius_Valleyproj/radiustable.csv', delimiter=',')
    masstable = np.loadtxt(r'C:/Users/madir/Radius_Valleyproj/masstable.csv', delimiter=',')
    
    return radtable, masstable
    
def Rcore_to_Mcore(Xiron, rcore):
    
    pos, coeff1, coeff2 = weights(Xiron)
    radtable, masstable = tables()
    
    
    f1 = interpolate.interp1d(masstable[:,pos-1], radtable[:,pos-1], bounds_error=False)
    f2 = interpolate.interp1d(masstable[:,pos], radtable[:,pos], bounds_error=False)
    massarr=np.logspace(-2,2,500)
    radtable1=f1(massarr)
    radtable2=f2(massarr)
    
    if coeff1 == 0: 
        radii = coeff2*radtable2
    else: 
        radii = coeff1*radtable1 + coeff2*radtable2
        
    #use linear interpolation to find the mass value for a given radius value as a function of Xiron
    massfunc = interpolate.interp1d(radii, massarr, bounds_error=False)
    mass = massfunc(rcore)
    return mass
    

def Mcore_to_Rcore(Xiron, mcore):
    
    pos, coeff1, coeff2 = weights(Xiron)
    radtable, masstable = tables()
    
    f1 = interpolate.interp1d(radtable[:,pos-1], masstable[:,pos-1], bounds_error=False)
    f2 = interpolate.interp1d(radtable[:,pos], masstable[:,pos], bounds_error=False)
    radarr=np.linspace(.1,10,500)
    masstable1 = f1(radarr)
    masstable2 = f2(radarr)
    
    if coeff1 == 0: 
        masses = coeff2*masstable2
    else: 
        masses = coeff1*masstable1 + coeff2*masstable2

    
    #use linear interpolation to find the radius value for a given mass value as a function of Xiron 
    radiusfunc = interpolate.interp1d(masses, radarr, bounds_error=False)
    radius = radiusfunc(mcore)
    
    return radius
    
def evaluate_X(Delta_Rrcb, Teq, Mcore, Tkh_Myr, Xiron):
    #solve for mass envelope frac using the radius of the radiative convect boundary
    
    #for density of core should I use planet.radius or core mass to core radius function? 
    rho_core=planet
    




    
    

    
    

    
    