#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 14:53:00 2018

@author: paul
"""

import numpy as np

class experiment:
    """
    object that defines one experiment
    """
    def __init__(self,duration,numff,flfrq,width_flfrq,fldur,rint,coupling_strength,moveparameter,vel):
        """
        initialize the experiment
        """
        self.flfrq     = flfrq 
        self.width_flfrq = width_flfrq
        self.fldur = fldur 
        
        if fldur>(0.5/flfrq):
            raise ValueError("flash duration can be maximum  half of the inverse flash frequency")
        
        self.rint = rint 
        self.coupling_strength = coupling_strength            
        self.vel =  vel
        self.number_fireflies    = numff
        self.duration =  duration
        self.dt = min(self.fldur/15/self.coupling_strength,self.fldur/15.) 
        self.iters = int(duration/self.dt)
        self.moveparameter = moveparameter
        self.flash_counter = np.zeros(self.iters)
        self.max_sim_flashs = np.zeros(self.iters)
        self.r_t = np.zeros(self.iters)
        
        # init fireflies
        self.firefly_list = [] 
        for i in range(self.number_fireflies):
            self.firefly_list.append(firefly(self.flfrq,self.width_flfrq,self.fldur,self.rint,self.coupling_strength,self.moveparameter,self.vel,self.dt))
    
    
    def run(self):
        """
        run the experiment for the chosen setting
        """
        t=0
        while True:
            self.print_status(t)
            self.time_stepper()
            flash_counter_t = self.count_flashs()
            t+=1
            
            # stop the integration if maximum number of iteration was reached or
            # if all fireflies are flashing
            if t>self.iters-1 or max(self.flash_counter)==self.number_fireflies:
                self.max_sim_flashs[t:]=self.number_fireflies
                self.r_t[t:] = 1.
                break
            
            self.flash_counter[t] = flash_counter_t
            self.max_sim_flashs[t] = max(self.flash_counter)
            self.r_t[t] =r_t_i(self.number_fireflies,self.firefly_list)
        
    
    def print_status(self,t):
        """ prints current sync state of the model"""
    
        test_print = t % 1000
        if test_print == 0.:
            print "Time:" + str(t*self.dt)+ "s, Sync:"+ str(self.max_sim_flashs[t]/self.number_fireflies *100.)+"%"
            
    def time_stepper(self):
        """performs one time step in integration for the model"""
        for i in range(self.number_fireflies):
            # time integration
            self.firefly_list[i].timeIntegration()
            
            # let the fireflies interact
            for j in range(self.number_fireflies):
                self.firefly_list[i].interact(self.firefly_list[j])
                
    def count_flashs(self):
        """ counts the number of flashing fireflies at the current time step"""
        flash_counter =0.
        for i in range(self.number_fireflies):
            flash_counter = flash_counter +  self.firefly_list[i].flash
        
        return flash_counter
        
         
            
class firefly():
    """
    object that defines the properties of one firefly
    """
    def __init__(self,flfrq,width_flfrq,fldur,rint,coupling_strength,moveparameter,vel,dt):
        """
        initialize the firefly properties
        """
        # initialize variables from settings
        self.flash_freq     = np.random.uniform(low=flfrq-width_flfrq,high=flfrq+width_flfrq)  
        self.flash_duration = fldur 
        self.thresold_flash = np.sin(np.pi/2.-self.flash_duration*self.flash_freq*np.pi) 
        self.r_interaction  = rint 
        self.coupling_strength = coupling_strength   
        self.moveparameter = moveparameter
        self.maxabsvel = vel
        self.maxvel = vel
        self.l = 1.0 
        self.dt = dt
        
        # random initialization of other variables
        self.x = np.random.uniform(low=0.,high=self.l)
        self.y = np.random.uniform(low=0.,high=self.l)
        self.intrinsic_omega    = 2*np.pi * self.flash_freq # np.random.uniform(low=0.9,high=1.0)
        self.omega    = self.intrinsic_omega
        self.theta    = np.random.uniform(low=0.,high=2*np.pi)
        self.flash    = flash(self.theta,self.thresold_flash)
        self.uvel     = np.random.uniform(low=-self.maxvel,high=self.maxvel)
        self.vvel     = np.random.uniform(low=-self.maxvel,high=self.maxvel)
    
    def timeIntegration(self):
        """
        integrate the model:
            - move the clock of a firefly
            - move the firefly on the grid (if moveparameter == True)
        """
        self.theta = self.theta + self.dt*self.omega % (2*np.pi)
        self.flash = flash(self.theta,self.thresold_flash)
        if self.moveparameter:        
            self.move()
            self.change_velocity()
        
    def interact(self,other_firefly):
        """
        interaction of a firefly with a other fireflies in its vicinity
        """
        Dx = other_firefly.x - self.x
        Dy =  other_firefly.y - self.y
        r = distance(Dx,Dy)
        
        if r!=0. and r<=self.r_interaction and self.flash==0. and other_firefly.flash==1.:
            self.theta = self.theta + np.sin(other_firefly.theta - self.theta) *self.coupling_strength*self.dt #% (2*np.pi) --> does not work caouse other.theta will also interact again!!! and if modlulo
            #other.theta +=self.K*self.dt
    
    def move(self):
        """
        move a firefly
        """
        dx = self.uvel*self.dt
        dy = self.vvel*self.dt    
        self.x = (self.x +dx) % self.l
        self.y = (self.y +dy) % self.l
        
    def change_velocity(self):
        """
        generate a new velocity vector
        """
        self.uvel     = self.uvel+np.random.uniform(low=-self.maxvel/10.,high=self.maxvel/10.)
        self.vvel     = self.vvel+np.random.uniform(low=-self.maxvel/10.,high=self.maxvel/10.)

# =============================================================================
# Some functions
# =============================================================================
# if a fireflie flashs than return 1.0, if not return 0.0
def flash(theta,thresold_flash):
    if np.sin(theta)>thresold_flash:
        return 1.
    else:
        return 0.
    
#calculate the distance between two points
def distance(x,y):
    return np.sqrt(x**2+y**2)

def r_t_i(nff,Fireflies):
    dummy = 0.
    for i in range(nff):
        dummy =  dummy + np.exp(1.0j * Fireflies[i].theta)
    return 1/float(nff) * np.abs(dummy)

def r_mean(r,T,T_0,iters):
    i_0 = int(T_0 * float(iters)/(T+T_0))
    return np.mean(r[i_0:])