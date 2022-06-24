# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#planet structure
import numpy as np
import constants as cs 
import planetstructure as ps 

#create objects to hold planet data



class Star():
    
    def __init__(self, mass, age, radius, Teff, mass_err=0, age_err=0, radius_err=0, Teff_err=0):
        self.mass = mass
        self.mass_err = mass_err
        self.age = age
        self.age_err = age_err
        self.radius = radius
        self.radius_err = radius_err
        self.Teff = Teff
        self.Teff_err = Teff_err

class Planet():
    
    def __init__(self,radius, period, isrocky, radius_err=0, period_err=0, Xiron=1/3, albedo=0):
        self.radius = radius
        self.radius_err = radius_err
        self.period = period
        self.period_err = period_err
        self.Xiron = Xiron
        self.albedo=albedo
        self.isrocky= isrocky
        self.calcMcore()
        
       
    def calcsemimajor(self, star):
            self.a = ((( cs.G *cs.Msun(star.mass)*(self.period*cs.d2s)**2)/(4*(np.pi)**2))**(1/3))/cs.Au2m
            
    def calcplanettemp(self, star):
            self.Teq = (1-self.albedo)**(1/4)*star.Teff*np.sqrt(cs.Rsun(star.radius)/(2*self.a*cs.Au2m))
            
    def calcMcore(self):
        if self.isrocky: 
            self.Mcore = ps.Rcore_to_Mcore(self)
            
        
            
    
class PlanetarySystem(): 
    
    def __init__(self, starparams, planetRockyparams, planetEnvparams, Mstar_err=0, Star_age_err=0, Rstar_err=0, Teff_err=0, Rrocky_err=0, Procky_err=0, Xironrocky=1/3, albedo_rocky=0, Renv_err=0, Penv_err=0, Xironenv=1/3, albedo_env=0 ):
        
        stardict={'mass_err':Mstar_err, 'age_err':Star_age_err, 'radius_err':Rstar_err, 'Teff_err':Teff_err}
        self.star=Star(*starparams, **stardict)
        
        planetRockydict={'radius_err':Rrocky_err, 'period_err':Procky_err, 'Xiron':Xironrocky, 'albedo':albedo_rocky}
        self.planetRocky=Planet(*planetRockyparams, isrocky=True, **planetRockydict)
        self.planetRocky.calcsemimajor(self.star)
        self.planetRocky.calcplanettemp(self.star)
        
        planetEnvdict={'radius_err':Renv_err, 'period_err':Penv_err, 'Xiron':Xironenv, 'albedo':albedo_env}
        self.planetEnv=Planet(*planetEnvparams, isrocky=False, **planetEnvdict)
        self.planetEnv.calcsemimajor(self.star)
        self.planetEnv.calcplanettemp(self.star)
        
        
        
    
            




            
            
            
            
            
            
            
    