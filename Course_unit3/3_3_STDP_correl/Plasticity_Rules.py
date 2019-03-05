# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 12:58:12 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg
"""

import numpy as np

class NN_STDP:
    
    ''''' this class defines an STDP object that allows the comparison between pre- and postsynaptic spikes
    and a matching weight update, based on the exponential STDP rule by Bi and Poo, 1998.
    The STDP is implemented in nearest-neighbour (NN) fashion.
    
    STDP_pre_post takes care of the 'pre then post' part of the window. 
    In the classical STDP this leads to LTP (when A_prepost is positive).
    
    STDP_post_pre takes care of the 'post then pre' part of the window. 
    In the classical STDP this leads to LTD (when A_postpre is negative).
    
    No change in the weight takes place if the difference in spike timing is zero.
    '''
    
    def STDP_pre_post(self,post_spike,pre_spikes_all,w_e_vec_tt,tau_prepost,A_prepost,w_max):

        # compare all most recent input spikes with this new post-spike.
        for syn in range(0,len(pre_spikes_all)):
            
            # get all spike-time differences for this synapse compared to the new post-spike
            diff = np.tile(post_spike,len(pre_spikes_all[syn])) - pre_spikes_all[syn] 
            diff_prepost = diff[np.where(diff>=0)] # take only spikes that happened before post.
        
            if np.size(diff_prepost > 0):
                Delta_t = - min(diff_prepost) # get nearest spike BEFORE post spike
                
                if Delta_t<0:
                    Delta_w_e=A_prepost*np.exp(Delta_t/tau_prepost)
                    #print('LTP happened')
                elif Delta_t>0:
                    print('Error in STDP code! Delta_t not consistent')
                    break
                else:
                    #print("pre_spike time = post_spike time")
                    Delta_w_e=0
                     
                # update the weight vector for this synapse, clip weight         
                w_e_vec_tt[syn]=w_e_vec_tt[syn]+Delta_w_e
                if w_e_vec_tt[syn]>w_max: # upper bound of the weight
                    w_e_vec_tt[syn]=w_max
                elif w_e_vec_tt[syn]<0: # lower bound of the weight
                    w_e_vec_tt[syn]=0
                         
            else:
                # no spikes from this synapse arrived before this post-spike.
                pass        

        return w_e_vec_tt
        
    def STDP_post_pre(self,post_spike,pre_spikes_all,tt,time_step_sim,w_e_vec_tt,tau_postpre,A_postpre,w_max):

        # compare all new pre-spikes with the latest post-spike.
        for syn in range(len(pre_spikes_all)):
            newest_spike = -1
            
            for sp in range(len(pre_spikes_all[syn])):
                # this synapse has at least one spike in this timestep:
                if pre_spikes_all[syn][sp]>=tt-time_step_sim and pre_spikes_all[syn][sp]<tt+time_step_sim:
                    newest_spike = pre_spikes_all[syn][sp] 
                    # choose only the spike closest to the post-spike (NN).
                    # ignore the subsequent spikes, even if they happen in same time bin and are more recent.
                    break 
                else: 
                    continue    
            
            if newest_spike >=0:
                # this synapse had at least one spike in this interval. (other synapses are ignored)
                    
                # get Delta t for its newest spike compared to last post-spike, 
                # which is the nearest spike AFTER the last post-spike
                Delta_t = newest_spike - post_spike  # should be positive
                
                if Delta_t>0:
                    Delta_w_e=A_postpre*np.exp(-Delta_t/tau_postpre)
                    #print('LTD happened')
                elif Delta_t<0:
                    print('Error in STDP code! Delta_t not consistent')
                    break
                else:
                    #print("pre_spike time = post_spike time")
                    Delta_w_e=0
                         
                # update the weight vector for this synapse, clip weights         
                w_e_vec_tt[syn]=w_e_vec_tt[syn]+Delta_w_e
                if w_e_vec_tt[syn]>w_max: # upper bound of the weight
                    w_e_vec_tt[syn]=w_max
                elif w_e_vec_tt[syn]<0: # lower bound of the weight
                    w_e_vec_tt[syn]=0
                                                
        return w_e_vec_tt    
        
            
        