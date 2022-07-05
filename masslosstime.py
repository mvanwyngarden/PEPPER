# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:29:40 2022

@author: madir
"""

import numpy as np 
import planetstructure as ps
import constants as cs
from scipy.optimize import minimize_scalar
import pdb


def calc_efficiency(Mcore, Rplanet): 
    '''Rplanet is in cm'''
    vesc = np.sqrt(2*cs.G*Mcore/Rplanet)
    #using varaible efficiency model from Owen & Wu 2017
    eta = 0.1*(vesc/15e3)**(-2)
    
    return eta

def calc_masslosstime_rocky(system, tmdot_Myr):
    
    #set constraints on envelope width for optimization function
    DR_rcb_min = 0.1 * cs.Rearth2cm(system.planetRocky.radius)
    DR_rcb_max = 5 * cs.Rearth2cm(system.planetRocky.radius)
    
    #optimization function that finds the value of DR_rcb that maximize mass loss time scale
    result=minimize_scalar(max_tmdot_objective, bounds=(DR_rcb_min, DR_rcb_max), args=(system, tmdot_Myr), method='bounded')
    
    if (result.success):
        DR_rcbmaximized = result.x
        #find maximized envelope mass fraction 
        X_max = ps.calc_X_adiabatic(cs.Rearth2cm(system.planetRocky.radius), DR_rcbmaximized, tmdot_Myr, system.planetRocky.Teq, system.planetRocky.Xiron)
        #find eta using Rp
        _, Rplanet = ps.calc_Rplanet(cs.Rearth2cm(system.planetRocky.radius), DR_rcbmaximized, tmdot_Myr, system.planetRocky.Teq, system.planetRocky.Xiron)
        eta = calc_efficiency(cs.Mearth2g(system.planetRocky.Mcore), Rplanet)
        
        masslosstime_max_rocky=masslosstime_eq(system, X_max, Rplanet, eta)
        
        system.planetRocky.X_max = X_max 
        system.planetRocky.DR_rcbmax = DR_rcbmaximized
        system.planetRocky.scaled_masslosstime_max = masslosstime_max_rocky
        
        return DR_rcbmaximized, masslosstime_max_rocky
    
    else: 
        print('Failed to find solution for maximum mass-loss timescale for rocky planet')
        return -1
        
def max_tmdot_objective(DR_rcb, system, Tkh_Myr): 
    #the mass loss timescale we want to maximize 
    
    
    X = ps.calc_X_adiabatic(cs.Rearth2cm(system.planetRocky.radius), DR_rcb, Tkh_Myr, system.planetRocky.Teq, system.planetRocky.Xiron)
    _, Rplanet = ps.calc_Rplanet(cs.Rearth2cm(system.planetRocky.radius), DR_rcb, Tkh_Myr, system.planetRocky.Teq, system.planetRocky.Xiron)

    eta = calc_efficiency(cs.Mearth2g(system.planetRocky.Mcore), Rplanet)
    func_to_max = X/ (eta*Rplanet**3)
    
    return 1/func_to_max    
    
    
def masslosstime_eq(system, X, Rplanet, eta):
    '''Rplanet is in cm'''
    
    masslosstime= X* system.planetRocky.Mcore**2 *system.planetRocky.a**2/ (Rplanet**3 *eta)
    return masslosstime



   
    
    
    
    
    
    