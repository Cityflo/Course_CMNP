# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:05:42 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

This class contains regularly spaced inputs (periodic input spikes)
for the LIF model.
"""

from __future__ import division
import random
import numpy as np

class Input_Regular: 
    
    #############################################
    # function for obtaining spike input,
    # with spike times given by argument 'spikes'.
    #############################################
    
    # g_tt_vec contains all the conductances from the last timestep
    # syn_func contains the function of the the conductances, of the synapse type (needed for Euler Integration)
    # weight is the synaptic strength
    # counter counts the number of inputs

    def inputs(self,g_tt_vec,time_step_sim,spikes,spiketrue,numb_syn,weight,syn_func,tau_syn,tt,counter,delta_t):
        
        from Euler_Method import Euler
        
        for syn in range(1,numb_syn+1):  # loop over synapses 
            
            # get spike times:                 
            if np.sum(spiketrue[np.where(np.tile(tt,len(spikes))==spikes)])>0: 
               
                # input is delivered   
                counter += 1
                g_tt_vec[syn-1]=g_tt_vec[syn-1]+weight[syn-1]
            
            # call function from Euler_Method to integrate the conductance
            g_1=Euler()
            g_tt_vec[syn-1]=g_1.euler_integration(syn_func, tau_syn, g_tt_vec[syn-1], tt, time_step_sim, delta_t) # gives us the value of the differential equation at certain time
        
        return [g_tt_vec,counter]

    