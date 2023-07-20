# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 10:51:59 2022

@author: madir
"""
import constants as cs
import pre_MS_luminosity as lumin
import numpy as np
import scipy as sc 
import planetstructure as ps
from scipy.optimize import brentq
import pdb



def calc_hill_radius(a, Mcore, star_mass):
    '''Returns the planet's hill radius in cm
    
    Parameters: 
        a - semimajor axis in cm
        Mcore - planet's core mass in Earth masses
        star_mass - star mass in Solar masses'''
    
    R_hill=a*(cs.Mearth2g(Mcore)/cs.Msun2g(3*star_mass))**(1/3)
    return R_hill
    
def calc_bondi_radius(Mcore, Tdisk):
    '''Returns the Bondi radius in cm
    
    Parameters: 
        Mcore - planet's core mass in Earth masses
        Tdisk - local disk temperature in K'''
    
    R_bondi = (cs.G*cs.Mearth2g(Mcore)*cs.mu/(cs.kb*Tdisk))
    
    return R_bondi

def calc_Rout(Mcore, a, star_mass, Tdisk): 
    '''Returns the bound radius of the envelope in cm
    
    Parameters: 
        Mcore - planet mass in Earth masses
        a - semimajor axis in cm
        star_mass - star mass in solar masses
        Tdisk - local disk temperature in K'''
    
    R_hill = calc_hill_radius(a, Mcore, star_mass)
    R_bondi = calc_bondi_radius(Mcore, Tdisk)
    
    #the bound radius is calculated by a factor f times the smaller of the Hill Radius or Bondi radius
    #f is a numerical factor set to 1.3 for gas-poor and dust-free atmospheres (Lee & Chiang, 2015)
    R_out=cs.f*min(R_hill, R_bondi)
    
    return R_out

def calc_disk_temp(star_mass, a, age):
    '''Returns the local disk temperature in K
    
    Parameters: 
        star_rad - Star radius in units of solar radii
        Teff - star's effective temperature in K
        a - semimajor axis in cm'''
    
    L = lumin.luminosity(star_mass, age)*3.839e33
    
    Tdisk = 550* (a/cs.Au2cm)**(-2/5)* (L/5.6e33)**(1/5)
    
    return Tdisk
    
def iso_atmos_mass_limit_eq(rho_disk, star_rad, Teff, Mcore, a, Rout, r, star_mass, age):
    '''The equation used to calculate the max envelope mass a planet can accrete given by its isothermal limit
    
    Parameters: 
        rho_disk - volumetric disk gas density
        star_rad - star radius in units of solar radii
        Teff - star's effective temperature in K
        Mcore - planet's core mass in Earth masses
        a - semimajor axis in cm
        Rout - planet radius in cm
        r - unit of integration, set to the integral bounds, either Rout or Rcore in cm
        star_mass - mass of the system's host star in solar masses
        age- age of the system in yrs'''
    
    Mcore=cs.Mearth2g(Mcore)
    
    Tdisk = calc_disk_temp(star_mass, a, age)
    c_atmos_sq = cs.kb*Tdisk/ cs.mu
    
    expi=sc.special.expi(cs.G*Mcore/(c_atmos_sq*r))
    
    p1= np.exp(-(cs.G*Mcore)/(c_atmos_sq*Rout))
    p2=c_atmos_sq*r*np.exp(cs.G*Mcore/(c_atmos_sq*r))
    p3=2*c_atmos_sq**2*r**2 + c_atmos_sq*cs.G*Mcore*r + cs.G**2*Mcore**2
    p4=cs.G**3 * Mcore**3 *expi
    
    M_iso_eq = 4*np.pi*rho_disk*(p1*(p2*p3-p4)/(6*c_atmos_sq**3))  
    
    return M_iso_eq 

def calc_X_iso(Rcore, star_rad, Teff, Mcore, a, star_mass, age): 
    '''Returns the maximum envelope mass fraction
    
    Parameters: 
        Rcore - planet radius in Earth radii
        star_rad - star radius in units of solar radii
        Teff - star's effective temperature in K
        Mcore - planet's core mass in Earth masses
        a - semimajor axis in cm
        star_mass - star mass in units of solar masses
        '''
    Tdisk = calc_disk_temp(star_mass, a, age)
    R_out = calc_Rout(Mcore, a, star_mass, Tdisk)
    M_iso = iso_atmos_mass_limit_eq(1, star_rad, Teff, Mcore, a, R_out, R_out,star_mass, age)-iso_atmos_mass_limit_eq(1, star_rad, Teff, Mcore, a, R_out, cs.Rearth2cm(Rcore), star_mass, age)
    
    return M_iso/cs.Mearth2g(Mcore)

def calc_min_mass_env(system, Rcore_env, Xiron, star_rad, Teff, a_env, star_mass, X_rocky, age):
    '''Returns the minimum core mass of the enveloped planet to be consistent with gas poor formation in Earth masses
    
    Parameters: 
        system - Planetary System class object
        Rcore_env - core radius of the enveloped planet in Earth masses
        Xiron - iron mass fraction
        star_rad - star radius in units of solar radii
        Teff - star's effective temperature in K
        a_env - semimajor axis of the enveloped planet in cm
        star_mass - star mass in units of solar masses'''
    
    if (Rcore_env>np.max(ps.table_interpolation(Xiron)[0])):
        Rmax = np.max(ps.table_interpolation(Xiron)[0])
        
    else: 
        Rmax = Rcore_env/1.1
        
    #set this to largest value in the grid
    Mcore_max = ps.Rcore_to_Mcore(Rmax, Xiron)
    
    Mcore_min = 0.1
    #pdb.set_trace()
    lg_min_mass_env=brentq(X_compare, np.log10(Mcore_min), np.log10(Mcore_max), args=(system, Rcore_env, star_rad, Teff, a_env, star_mass, X_rocky, age))
    min_mass_env= 10**lg_min_mass_env
    
    return min_mass_env
    
def X_compare(lg_Mcore, system, Rcore_env, star_rad, Teff, a_env, star_mass, X_rocky, age):
    '''Returns the difference between the envelope mass fraction of the rocky planet and the envelope mass fraction of the enveloped planet
    
    Parameters: 
        Mcore - core mass of the enveloped planet 
        system - Planetary System class object
        Rcore_env - core radius of the enveloped planet in Earth masses
        star_rad - star radius in units of solar radii
        Teff - star's effective temperature in K
        a_env - semimajor axis of the enveloped planet in cm
        star_mass - star mass in units of solar masses
        X_rocky - envelope mass fraction of rocky planet'''
    
    Mcore = 10**lg_Mcore
    
    X_env= calc_X_iso(Rcore_env, star_rad, Teff, Mcore, a_env, star_mass, age)
    
    return X_rocky-X_env
    
    
    


    
