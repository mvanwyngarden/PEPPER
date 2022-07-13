# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#planet structure
import numpy as np
import constants as cs 
import planetstructure as ps 
import masslosstime as ml
import masslosstime_env as ml_env
import pdb

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
    
    def __init__(self, radius, period, isrocky, radius_err=0, period_err=0,  Xiron=1/3, albedo=0):
        self.radius = radius   #make numpy array with vals from a guassian distribution whose standard dev is rad err
        self.radius_err = radius_err
        self.period = period
        self.period_err = period_err
        self.isrocky= isrocky
        self.Xiron = Xiron
        self.albedo=albedo
        
        
        
    
    
    
class PlanetarySystem(): 
    
    
    def __init__(self, starparams, planetRockyparams, planetEnvparams,  N=1000, Tkh_Myr=100, Mstar_err=0, Star_age_err=0, Rstar_err=0, Teff_err=0, Rrocky_err=0, Procky_err=0, Xironrocky=1/3, albedo_rocky=0, Renv_err=0, Penv_err=0, Xironenv=1/3, albedo_env=0 ):
        
        stardict={'mass_err':Mstar_err, 'age_err':Star_age_err, 'radius_err':Rstar_err, 'Teff_err':Teff_err}
        self.star=Star(*starparams, **stardict)
        
        planetRockydict={'radius_err':Rrocky_err, 'period_err':Procky_err, 'Xiron':Xironrocky, 'albedo':albedo_rocky}
        self.planetRocky=Planet(*planetRockyparams, isrocky=True, **planetRockydict)
        
       
        planetEnvdict={'radius_err':Renv_err, 'period_err':Penv_err, 'Xiron':Xironenv, 'albedo':albedo_env}
        self.planetEnv=Planet(*planetEnvparams, isrocky=False, **planetEnvdict)
        
        errs=np.array([Mstar_err, Star_age_err, Rstar_err, Teff_err, Rrocky_err, Procky_err, Renv_err, Penv_err])
        self.N =int(N) if np.any(errs>0) else 1
        
        self.Tkh_Myr = Tkh_Myr
     
        self.sample_params()
        self.calc_min_mass_env()
        
    def calcsemimajor(self,star_mass, planet_period):
        a = ((( cs.G *cs.Msun2g(star_mass)*(planet_period*cs.d2s)**2)/(4*(np.pi)**2))**(1/3))
        return a
            
    def calcplanettemp(self,albedo, star_Teff, star_radius, star_mass, planetRocky_period):
        a = ((( cs.G *cs.Msun2g(star_mass)*(planetRocky_period*cs.d2s)**2)/(4*(np.pi)**2))**(1/3))
        Teq = (1-albedo)**(1/4)*star_Teff*np.sqrt(cs.Rsun2cm(star_radius)/(2*a))
        return Teq
            

       
    def sample_params(self):
        self.planetRocky.Rcore_samp = np.random.normal(self.planetRocky.radius, self.planetRocky.radius_err, self.N)
        self.planetEnv.Rcore_samp = np.random.normal(self.planetEnv.radius, self.planetEnv.radius_err, self.N)
        self.star.radius_samp = np.random.normal(self.star.radius, self.star.radius_err, self.N)
             
        self.planetRocky.period_samp = np.random.normal(self.planetRocky.period, self.planetRocky.period_err, self.N)
        self.planetEnv.period_samp = np.random.normal(self.planetEnv.period, self.planetEnv.period_err, self.N)
             
        self.star.mass_samp = np.random.normal(self.star.mass, self.star.mass_err, self.N)
    
        self.star.age_samp = np.random.normal(self.star.age, self.star.age_err, self.N)
    
        self.star.Teff_samp = np.random.normal(self.star.Teff, self.star.Teff_err, self.N)
         
        
    def calc_min_mass_env(self):
        
        self.planetEnv.minMcore_samps = np.zeros(self.N)
        for i in range (self.N):
            
            
            self.planetRocky.a = self.calcsemimajor(self.star.mass_samp[i], self.planetRocky.period_samp[i])
            self.planetEnv.a = self.calcsemimajor(self.star.mass_samp[i], self.planetEnv.period_samp[i])
            
            self.planetRocky.Teq = self.calcplanettemp(self.planetRocky.albedo, self.star.Teff_samp[i], self.star.radius_samp[i], self.star.mass_samp[i], self.planetRocky.period_samp[i])
            self.planetEnv.Teq = self.calcplanettemp(self.planetEnv.albedo, self.star.Teff_samp[i], self.star.radius_samp[i], self.star.mass_samp[i], self.planetEnv.period_samp[i])
            
            self.planetRocky.Mcore=ps.Rcore_to_Mcore(self.planetRocky.Rcore_samp[i], self.planetRocky.Xiron)
            try: 
                
                ml.calc_masslosstime_rocky(self, self.planetRocky.Rcore_samp[i], self.planetRocky.Mcore, self.planetRocky.a, self.planetRocky.Teq, self.planetRocky.Xiron, self.Tkh_Myr)
            except: 
                pass
            try:
                
                min_mass = ml_env.min_mass_env(self.planetEnv.Rcore_samp[i], self.planetEnv.a, self.planetEnv.Teq, self.star.age_samp[i], self.Tkh_Myr, self.planetEnv.Xiron, self.planetRocky.scaled_masslosstime_max)
            except: 
                pass
            
            self.planetEnv.minMcore_samps[i]= min_mass
        
    
            




            
            
            
            
            
            
            
    