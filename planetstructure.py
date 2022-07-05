# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
from scipy.integrate import quad
from bisect import bisect_left
import constants as cs
from scipy.optimize import fsolve
import pdb

def weights(Xiron):
    '''Returns the coefficients needed to find the weighted average of two mass/radius columns
       
       Parameters: Xiron'''
    
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
    
def Rcore_to_Mcore(Rcore, Xiron):
    '''Returns the mass of a planet given a radius value in Earth radii. 
    Mass-radius relationships derived from Li Zeng's models of Fe-MgSiO3 planets'''
    
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
    mass = massfunc(Rcore)
    return mass
    

def Mcore_to_Rcore(Mcore, Xiron):
    '''Returns the radius of a planet given an Xiron fraction and a mass value in Earth masses. 
    Mass-radius relationships derived from Li Zeng's models of Fe-MgSiO3 planets'''
    
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
    radius = radiusfunc(Mcore)
    
    return radius
    
##Envelope structure 

def integrand1(x):

    return x*(1/x-1)**(1/(cs.gamma-1))

def integrand2(x):

    return x**2*(1/x-1)**(1/(cs.gamma-1))

#integrates eq 5 in Owen & Wu from Rc/Rrcb to 1 

def get_I1(Rcore, Rrcb):
    
 
    I1 = quad(integrand1, Rcore/Rrcb, 1)[0]
    return I1

def get_I2(Rcore, Rrcb):
    #make input args same units
   
    #print (Rcore/Rrcb, Rcore, Rrcb)
    I2 = quad(integrand2, Rcore/Rrcb, 1)[0]
   
    
    return I2 

def get_I2_I1(Rcore, Rrcb):
    
    
    ratio = get_I2( Rcore, Rrcb)/get_I1( Rcore, Rrcb)
    return ratio 

def calc_X_adiabatic(Rcore, DR_rcb, Tkh_Myr, Teq, Xiron): 
    ''''
    Calculates the envelope mass fraction X
    
    Parameters: 
        Rcore - cm
        DR_rcb - the width of the adiabatic portion of the atmosphere in cm (float) 
        Tkh_Myr - cooling timescale in Myr
        Teq- planet equilibrium temp in K
        Xiron - '''
    
   
      
    Rrcb = DR_rcb + Rcore
   
    #sound speed of atmosphere
    c_atmos_sq = cs.kb * Teq / cs.mu
    
    adiabatic_grad = (cs.gamma-1)/cs.gamma
    
    #average core density 
    rho_core = calc_coredensity(Rcore, Xiron)
    
    #using eq 13 from Owen & Wu 2017, factor out X
    
    rho_rcb_factoroutx = calc_rho_rcb( Rcore, np.log10(DR_rcb), 1, Tkh_Myr, Teq, Xiron)
    
    
    #X=Menv/Mcore
    #by factoring out X from eq 13 we can replace the planet density in eq 4 with rho_rcb times x**(1+alpha)
    #using eq 4 from Owen & Wu, divide out Mcore, resulting in an expression of Menv/Mcore 
    Menv_Mcore = 3 * (Rrcb/Rcore)**3 * (rho_rcb_factoroutx/rho_core) * (adiabatic_grad *\
                    (cs.G * cs.Mearth2g(Rcore_to_Mcore(cs.cm2Rearth(Rcore), Xiron))) / (c_atmos_sq * Rrcb))**(1/(cs.gamma-1)) * get_I2(Rcore, Rrcb)
   
    #because Menv/Mcore=X we can solve for X 
    X = Menv_Mcore**(1/(1+(1/(1+cs.alpha))))
   
    
    return X
    

def calc_coredensity(Rcore, Xiron): 
    '''Returns the average density of the core in cgs
    Rcore - cm'''
    
    rho_core_calc = cs.Mearth2g(Rcore_to_Mcore(cs.cm2Rearth(Rcore), Xiron))/ ( (4/3)*np.pi*Rcore**3)
   
    return rho_core_calc



def calc_rho_rcb(Rcore, lg_DR_rcb, X, Tkh_Myr, Teq, Xiron): 
    '''Calculate the density at the rcb
    
    Parameters: 
        PlanetarySystem object
        DR_rcb- the width of the adiabatic portion of the atmosphere in cm (float) 
        X - envelope mass fraction (float)
       Tkh_Myr - cooling timescale in Myr'''
    
    
    #height of the rcb
    Rrcb = 10**lg_DR_rcb + Rcore
    
    I2_I1 = get_I2_I1(Rcore, Rrcb)
    
    Tkh_sec= cs.Myr2sec(Tkh_Myr)
    
    
    rho_rcb = (cs.mu/cs.kb) * (I2_I1 * 64 * np.pi * cs.sigma * Teq**(3-cs.alpha-cs.beta) * Rrcb * Tkh_sec / 
                            (3 * cs.kappa0 * cs.Mearth2g(Rcore_to_Mcore(cs.cm2Rearth(Rcore), Xiron))* X))**(1/(1+cs.alpha))
    
    return rho_rcb
 
    
def calc_Rplanet(Rcore, DR_rcb, Tkh_Myr, Teq, Xiron):
     '''Returns the radius of the planet over the radius to the rcb and returns the radius of the planet 
     Parameters: 
         PlanetarySystem object
         DR_rcb- the width of the adiabatic portion of the atmosphere in cm (float)
         Tkh_Myr- cooling timescale in Myr'''
     
     #calculate the denisty at the photosphere
     #pressure at the rcb is used as an approx for photospheric pressure 
     #because the isothermal layer is assumed to be thin
     
     #height of the rcb
     
     Rrcb = DR_rcb + Rcore
     
     Mcore = Rcore_to_Mcore(cs.cm2Rearth(Rcore), Xiron)
     
     
     X=calc_X_adiabatic( Rcore, DR_rcb, Tkh_Myr, Teq, Xiron)
     
     #pressure at the photosphere
     
     #do I need to change whether you insert Rrcb here? 
     _, rho_phot, H = calc_photo_pressure( Rrcb, Mcore, Teq)
    
    
     #f=Rp/Rrcb
     rho_rcb = calc_rho_rcb(Rcore, np.log10(DR_rcb), X, Tkh_Myr, Teq, Xiron)
     Rp_Rrcb = 1 + (H/Rrcb) * np.log(rho_rcb/rho_phot)
     
    
     
     if Rp_Rrcb < 1: 
         print('Error: rho_rcb < rho_phot')
         
     else: 
        Rplanet = Rp_Rrcb * Rrcb
        return Rp_Rrcb, Rplanet
    
 
def calc_X_isothermal(Rp, Mcore, Teq, Xiron): 
    #we'll need this equation if the adiabatic layer is less than 1 scale height
    #figure out how/where you want to call this
    Rcore = cs.Rearth2cm(Mcore_to_Rcore(Mcore, Xiron))
    
    _, rho_phot, H = calc_photo_pressure(Rp, Mcore, Teq)
    
    #find mass in radiative layer
    rho_atmos = rho_phot * np.exp(Rp/H * (Rp/Rcore-1))
    
    Menv = rho_atmos * 4 * np.pi * Rcore**2 * H
    
    X = Menv / (cs.Mearth2g(Mcore))
    
    return X, Rp/Rcore, Rp
    
    
def calc_photo_pressure( R, Mcore, Teq): 
    #here you will either insert Rrcb or Rp depending on if you need the adiabatic or isothermal layer 
    pressure_phot = (2/3 * (cs.G*cs.Mearth2g(Mcore)/(R**2.*cs.kappa0*Teq**cs.beta)))**(1/(1+cs.alpha))
    
    rho_phot_calc = (cs.mu/cs.kb) * pressure_phot / Teq
    
    H = cs.kb * Teq * R** 2 / ( cs.mu * cs.G * cs.Mearth2g(Mcore))
    
    return pressure_phot, rho_phot_calc, H
   
    
def Rplanet_solver(Rp, Mcore, planet_age, Teq, Xiron):
    
    #this will only get called when evaluating the gaseous planet because we will use the 
    #current radius of the planet as Rp 
    #we will then make sure that the Rplanet returned converges with the current Rp input by the user 
    Rcore = Mcore_to_Rcore(Mcore, Xiron)
    
    if (Rp < Rcore):
        return -1 
    
    Rcore = cs.Rearth2cm(Rcore)
    
    DR_rcb_guess = cs.Rearth2cm(Rp) - Rcore
    lg_DR_rcb_guess = np.log10(DR_rcb_guess)
    
    #pdb.set_trace()
    #find the DR_rcb where the calculated Rplanet is equal to the observed Rplanet 
    lg_DR_rcb_sol = fsolve(Rplanet_solver_eq, lg_DR_rcb_guess, args=(Rp, Mcore, Teq, planet_age, Xiron))
    
    H = calc_photo_pressure( cs.Rearth2cm(Rp), Mcore, Teq)[2]
    
    DR_rcb_sol = 10**lg_DR_rcb_sol
  
    if DR_rcb_sol > H: 
        X = calc_X_adiabatic(Rcore, DR_rcb_sol, planet_age, Teq, Xiron)
        Rp_Rrcb, Rplanet = calc_Rplanet(Rcore, DR_rcb_sol, planet_age, Teq, Xiron)
        return X, Rp_Rrcb, Rplanet
    else: 
        X, Rp_Rcore, Rplanet = calc_X_isothermal( cs.Rearth2cm(Rp), Mcore, Teq, Xiron)
        return X, Rp_Rcore, Rplanet
  #this returns the X and Rplanet at the current age of the system  
    

def Rplanet_solver_eq(lg_DR_rcb, Rp, Mcore, Teq, planet_age, Xiron): 
    
    #evaluate the envelope mass fraction and atmosphere structure for this mass/radius guess
    Rcore = cs.Rearth2cm(Mcore_to_Rcore(Mcore, Xiron))
   
    DR_rcb = 10**lg_DR_rcb
    print(DR_rcb)
    #does Rplanet return a value in cm? 
    Rp_Rrcb, Rplanet = calc_Rplanet(Rcore, DR_rcb, planet_age, Teq, Xiron)
    
    return cs.Rearth2cm(Rp) - Rplanet






    
    
    
    
    
    
     
      

    
    
    
    




    
    

    
    

    
    