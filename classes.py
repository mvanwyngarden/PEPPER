# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#planet structure
import numpy as np
import constants as cs 

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
    
    def __init__(self,radius, period, star, radius_err=0, period_err=0, Xiron=1/3, albedo=0):
        self.radius = radius
        self.radius_err = radius_err
        self.period = period
        self.period_err = period_err
        self.Xiron = Xiron
        self.albedo=albedo
        self.semimajor(star)
        self.planettemp(star)
       
    def semimajor(self, star):
            self.a = ((( cs.G *star.mass*cs.Msun*(self.period*cs.d2s)**2)/(4*(np.pi)**2))**(1/3))/cs.Au2m
            
    def planettemp(self, star):
            self.Teq = (1-self.albedo)**(1/4)*star.Teff*np.sqrt(star.radius*cs.Rsun/(2*self.a*cs.Au2m))
            
    def isrocky(self): 
        if self.radius < 1.8: 
            #rocky planet
            return True 
        else: 
            #enveloped planet
            return False 
            
            
star= Star(1, 2, 1, 5700)
planet1=Planet(1, 365, star, albedo=1)




            
            
            
            
            
            
            
    