# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 15:13:13 2022

@author: madir
"""

import numpy as np
import planetstructure as ps 
import constants as cs 
import masslosstime as ms
from scipy.optimize import minimize_scalar
from scipy.optimize import brentq
import pdb

def mass_limits_env(system, Teq, planet_age, Tkh_Myr, Xiron): 
    '''Returns the minimum and maximum mass values for which a mass-loss timescale can be calculated.
    
    Parameters: 
        system- Planetary system containing a rocky and enveloped planet and a star
        Teq - enveloped planet equilibrium temp, K
        planet_age - age of the system, Myr
        Tkh_Myr - cooling timescale, Myr
        Xiron - iron mass fraction
        '''
    
    #we want to make sure a solution exists, so we will find the maximum mass-loss timescale for the enveloped
    #and make sure that it goes above the maximum timescale for the rocky planet 
    
    #find the maximum core mass to check up to 
    Rcore = system.planetEnv.radius / 1.1 
    Mcore_max = ps.Rcore_to_Mcore(Rcore, Xiron)
    
    Mcore_min = 0.1
    #check that a minimum timescale exists for this low of a core mass 
    #if not, increase the min core mass until a solution exists 
    
    while Mcore_min < Mcore_max: 
        
        solution = masslosstime_env_minimize(Mcore_min, system, Teq, planet_age, Tkh_Myr, Xiron)
        
        if solution < 0 :
            Mcore_min += 0.1
            
        else: 
            if (1/solution < system.planetRocky.scaled_masslosstime_max): 
                #this can be used as a lower mass bound
                break
            
            else: 
                return -1 
            
    if solution < 0: 
        #could not solve for a lower mass bound 
        #no mass exists for which the mass loss timescale is less than the maximized timescale for the rocky
        return -2
    
    masslosstime_env_min = 1/masslosstime_env_minimize(Mcore_min, system, Teq, planet_age, Tkh_Myr, Xiron)
    
    #now find the upper bound for the mass value 
    #this returns the mass value that maximizes the mass loss timescale
    result = minimize_scalar(masslosstime_env_minimize, bounds = (Mcore_min, Mcore_max), args=(system, Teq, planet_age, Tkh_Myr, Xiron), method='bounded')
    
    if result.success: 
        
        Mcore_max = result.x
        
        masslosstime_env_max = 1/(masslosstime_env_minimize(Mcore_max, system, Teq, planet_age, Tkh_Myr, Xiron))
        
        #if the max timescale for the env is less than the rocky there is no solution 
        #because the env would have already lost its envelope
        
        if (masslosstime_env_max < system.planetRocky.scaled_masslosstime_max): 
            return -4 
        
    else: 
        #could not find max Mcore for the enveloped planet
        return -5 
    
    Mcore_max = float(Mcore_max)
    
    return Mcore_min, Mcore_max

    
def min_mass_env(system, Teq, planet_age, Tkh_Myr, Xiron, Mcore_min, Mcore_max, masslosstime_want): 
    '''Returns the minimum mass for the enveloped planet, this occurs when the mass-loss timescale for the enveloped planet equals the maximized
    mass-loss timescale for the rocky planet
    
    Parameters: 
        system- Planetary system containing a rocky and enveloped planet and a star
        Teq - enveloped planet equilibrium temp, K
        planet_age - age of the system, Myr
        Tkh_Myr - cooling timescale, Myr
        Xiron - iron mass fraction
        Mcore_min - minimum possible mass returned by mass_limits_env
        Mcore_max - maximum possible mass returned by mass_limits_env
        masslosstime_want - the maximized mass-loss timescale for the rocky planet'''
        
        
    #now we use these mass limits as bounds for our root solver
    #need to find when timescale for gas = max timescale for rocky
    #use log values to make the solver behave better 
    
    lg_min_mass_sol = brentq(masslosstime_env, np.log10(Mcore_min), np.log10(Mcore_max), args=(system, Teq, planet_age, Tkh_Myr, Xiron, masslosstime_want))
    min_mass_sol = 10**lg_min_mass_sol
    
    return min_mass_sol
    

    
def calc_Rplanet_scaledtime_env(system, Mcore, Teq, planet_age, Tkh_Myr, Xiron): 
    '''Returns the enveloped planet radius and envelope mass fraction at the time we wish to scale to. This is usually set at the cooling timescale of 100 Myr
    
    Parameters: 
        system- Planetary system object containing a rocky and enveloped planet and a star
        Mcore - planet mass value in Earth masses
        Teq - enveloped planet equilibrium temp, K
        planet_age - age of the system, Myr
        Tkh_Myr - cooling timescale, Myr
        Xiron - iron mass fraction'''
        
    Rplanet_now = system.planetEnv.radius
    #find X, f, and Rplanet returned by Rp_solver using Rplanet_now
    
    X, Rp_Rc, Rp_solved = ps.Rplanet_solver(Rplanet_now, Mcore, planet_age, Teq, Xiron )
    #we should compare the Rp_solved with the Rplanet observed so that we know Rp solver converged 
   
    if (np.abs(Rp_solved - cs.Rearth2cm(Rplanet_now)) / cs.Rearth2cm(Rplanet_now) > 1e-4): 
        return -6 
    
    #we know that the solver returned the correct values so we can use this value of X going forward
    
    #find the radius of the planet at the time we wish to scale to 
    Rplanet_tscale, DR_rcb = ps.solve_envelope_structure(X, Mcore, Tkh_Myr, Teq, Xiron)
    
    #check that the solver found the correct envelope structure by analytically calculating X 
    
    Rcore = ps.Mcore_to_Rcore(Mcore, Xiron)
    
    X_analytic = ps.calc_X_adiabatic(Rcore, DR_rcb, Tkh_Myr, Teq, Xiron)
    
    if (np.abs(X_analytic - X)/X > 1e-4): 
        return -7 
    
    else: 
        return float(Rplanet_tscale), float(X)

   
def masslosstime_env(lg_Mcore, system, Teq, planet_age, Tkh_Myr, Xiron, masslosstime_want): 
    '''Returns the mass loss timescale for the enveloped planet minus a given mass loss timescale. 
    The mass loss timescale for the enveloped planet is calculated using the scaled radius value. 
    
    Parameters: 
        lg_Mcore - log planet mass in Earth masses
        system - Planetary system object containing a rocky and enveloped planet and a star
        Teq - enveloped planet equilibrium temp, K
        planet_age - age of the system, Myr
        Tkh_Myr - cooling timescale, Myr
        Xiron - iron mass fraction
        masslosstime_want - a mass loss timescale, eventually this is set to the maximized mass loss timescale for the rocky planet'''
     
    Mcore = 10**lg_Mcore
     
    Rplanet_tscale, X = calc_Rplanet_scaledtime_env(system, Mcore, Teq, planet_age, Tkh_Myr, Xiron)
     
    eta = ms.calc_efficiency(cs.Mearth2g(Mcore), Rplanet_tscale)
     
    masslosstime_env = ms.masslosstime_eq(Mcore, system.planetEnv.a, X, Rplanet_tscale, eta)
     
    return masslosstime_env - masslosstime_want
 
    
def masslosstime_env_minimize(Mcore, system, Teq, planet_age, Tkh_Myr, Xiron): 
    '''Returns the inverse of the mass loss timescale for the enveloped planet
    
    Parameters:
        Mcore - enveloped planet mass in Earth masses
        system - Planetary system object containing a rocky and enveloped planet and a star
        Teq - enveloped planet equilibrium temp, K
        planet_age - age of the system, Myr
        Tkh_Myr - cooling timescale, Myr
        Xiron - iron mass fraction'''
        
    masslosstime_want = 0
    #set to 0 because we aren't ready to compare the enveloped mass loss timescale to the timescale of the rocky planet 
    
    lg_Mcore = np.log10(Mcore)
    
    masslosstime_env_min = masslosstime_env(lg_Mcore, system, Teq, planet_age, Tkh_Myr, Xiron, masslosstime_want)
    
    return 1/masslosstime_env_min 








    
    