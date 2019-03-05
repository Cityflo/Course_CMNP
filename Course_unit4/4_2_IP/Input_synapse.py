# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 14:05:42 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg
"""

import random

class Synapse:
    
    ###############
    '''
    Class for calculating all the excitatory or inhibitory conductances for a given timestep.
    It can be called for either excitatory or inhibitory inputs.
    The spike trains for each input have been pre-calculated,
    and simply inserted via 'spike_trains_complete'.
    
    Requires Euler integration class.
    '''
    ###############
    
        
    def inputs_calc(self,g_tt_vec,time_step_sim,numb_syn,w_vec_tt,syn_func,tau_syn,tt,delta_t,spike_trains_complete):
        
        from Euler_Method import Euler # import Euler to perform the Euler integration method
        
        for syn in range(1,numb_syn+1): # loop through synapses
            for sp in range(1,len(spike_trains_complete[syn-1])): # loop through spikes
        
                if spike_trains_complete[syn-1][sp-1]>=tt-time_step_sim and spike_trains_complete[syn-1][sp-1]<tt+time_step_sim:
                    g_tt_vec[syn-1]=g_tt_vec[syn-1]+w_vec_tt[syn-1]

            # call function from Euler_Method to integrate the conductance
            g_1=Euler()
            g_tt_vec[syn-1]=g_1.euler_integration(syn_func, tau_syn, g_tt_vec[syn-1], tt, time_step_sim, delta_t) # gives us the value of the differential equation at certain time
        
        return [g_tt_vec]      
    