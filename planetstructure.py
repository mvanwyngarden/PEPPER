# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
from scipy.integrate import quad
from bisect import bisect_left
import constants as cs
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
    
def Rcore_to_Mcore(planet):
    
    pos, coeff1, coeff2 = weights(planet.Xiron)
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
    mass = massfunc(planet.radius)
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
    
##Envelope structure 

def integrand1(x, gamma):

    return x*(1./x-1.)**(1./(gamma-1.))

def integrand2(x, gamma):

    return x**2.*(1./x-1.)**(1./(gamma-1.))

#integrates eq 5 in Owen & Wu from Rc/Rrcb to 1 

def get_I1(system, Rrcb, gamma):
    
    Rcore = cs.Rearth2cm(system.planetRocky.radius)
    I1 = quad(integrand1, Rcore/Rrcb, 1, args=gamma)[0]
    return I1

def get_I2(system, Rrcb, gamma):
    
    Rcore = cs.Rearth2cm(system.planetRocky.radius)
    I2 = quad(integrand2, Rcore/Rrcb, 1, args=gamma)[0]
    return I2 

def get_I2_I1(system, Rrcb, gamma):
    
    
    ratio = get_I2(system, Rrcb, gamma)/get_I1(system, Rrcb, gamma)
    return ratio 

def calc_X(system, DR_rcb, Tkh_Myr, gamma): 
    ''''
    Calculates X given planet parameters'''
    
    Rcore = cs.Rearth2cm(system.planetRocky.radius) 
    Rrcb = DR_rcb + Rcore
    
    #sound speed of atmosphere
    c_atmos_sq = cs.kb * system.planetRocky.Teq / cs.mu
    
    adiabatic_grad = (gamma-1)/gamma
    
    #average core density 
    rho_core = calc_coredensity(system)
    
    #using eq 13 from Owen & Wu 2017, factor out X
    
    rho_rcb_factoroutx = calc_rho_rcb(system, DR_rcb, 1, Tkh_Myr)
    
    
    #X=Menv/Mcore
    #by factoring out X from eq 13 we can replace the planet density in eq 4 with rho_rcb times x**(1+alpha)
    #using eq 4 from Owen & Wu, divide out Mcore, resulting in an expression of Menv/Mcore 
    Menv_Mcore = 3 * (Rrcb/Rcore)**3 * (rho_rcb_factoroutx/rho_core) * (adiabatic_grad *\
                    (cs.G * cs.Mearth(system.planetRocky.Mcore)) / (c_atmos_sq * Rrcb))**(1/(gamma-1)) * get_I2(system, Rrcb, gamma)
   
    #because Menv/Mcore=X we can solve for X 
    pdb.set_trace()
    X = Menv_Mcore**(1/(1+(1/(1+cs.alpha))))
    
    return X
    

def calc_coredensity(system): 
    #find average density of the core 
    
    rho_core_calc = cs.Mearth(system.planetRocky.Mcore)/ ( (4/3)*np.pi*(cs.Rearth2cm(system.planetRocky.radius))**3)
   
    return rho_core_calc



def calc_rho_rcb(system, DR_rcb, X, Tkh_Myr): 
    '''calculate the density at the rcb using your planet system and DR_rcb in cm'''
    #within the function make all constants into cgs units
    #then convert rho back into SI
    
    Rcore= cs.Rearth2cm(system.planetRocky.radius)
    
    #height of the rcb
    #DR_rcb is in Earth radii
    Rrcb = DR_rcb + Rcore

    I2_I1 = get_I2_I1(system, Rrcb, cs.gamma)
    
    Tkh_sec= cs.Myr2sec(Tkh_Myr)
    
    
    rho_rcb = (cs.mu/cs.kb) * (I2_I1 * 64 * np.pi * cs.sigma * system.planetRocky.Teq**(3-cs.alpha-cs.beta) * Rrcb * Tkh_sec / 
                            (3 * cs.kappa0 * cs.Mearth(system.planetRocky.Mcore)* X))**(1/(1+cs.alpha))
    
    return rho_rcb
 
#ok so I believe that this eq is correct
#but calc_x is not
    
    
    
    




    
    

    
    

    
    