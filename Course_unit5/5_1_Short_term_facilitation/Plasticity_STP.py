# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 12:58:12 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg
"""

import numpy as np

class STPlasticity:
    
    ''''' This class defines a short-term plasticity (STP) object.
    It is a type of short-term change in the synaptic strength due to
    changes in the presynapse.
    Here, only short-term facilitation (STF) is included.
    
    '''
    
    def STP(self,pre_spikes_all,tt,time_step_sim,tau_f,tau_d,U,w_fixed,w_e_vec_tt,u_vec,x_vec):
           
        # did the synapse have a pre-spike in this interval?   
        for syn in range(len(pre_spikes_all)):
            
            spike_or_not=0    
            for sp in range(len(pre_spikes_all[syn])):
                # this synapse has at least one spike in this timestep:
                if pre_spikes_all[syn][sp]>=tt-time_step_sim and pre_spikes_all[syn][sp]<tt+time_step_sim:
                    newest_spike = pre_spikes_all[syn][sp] 
                    
                    # Short-term faciliation
                    u_vec[syn]=u_vec[syn]-u_vec[syn]/tau_f+U*(1-u_vec[syn])
            
                    # weight change
                    w_e_vec_tt[syn]=w_fixed*u_vec[syn]*x_vec[syn]
                    
                    spike_or_not=1
                    
                    break 
                else: 
                    continue
            
            if spike_or_not==0:
                
            # recovery of the STF variable for this synapse
                # Short-term faciliation
                u_vec[syn]=u_vec[syn]-u_vec[syn]/tau_f
                
        return [w_e_vec_tt,u_vec,x_vec]
        
            
        