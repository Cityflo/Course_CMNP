# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:33:59 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

This code runs a LIF-neuron model with excitatory and inhibitory synaptic inputs.
The synaptic inputs arrive periodically.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

from Euler_Method import Euler # import Euler to perform the Euler integration method
from Input_Regular import Input_Regular # import class where the excitatory and inhibitory input is calculated

import Parameters_Int_and_Fire

tau_mem       = Parameters_Int_and_Fire.tau_mem
E_leak        = Parameters_Int_and_Fire.E_leak
E_e           = Parameters_Int_and_Fire.E_e
E_i           = Parameters_Int_and_Fire.E_i
V_reset       = Parameters_Int_and_Fire.V_reset
V_thresh      = Parameters_Int_and_Fire.V_thresh
t_0           = Parameters_Int_and_Fire.t_0
t_max         = Parameters_Int_and_Fire.t_max
time_step_sim = Parameters_Int_and_Fire.time_step_sim
numb_exc_syn  = Parameters_Int_and_Fire.numb_exc_syn
numb_inh_syn  = Parameters_Int_and_Fire.numb_inh_syn
tau_e         = Parameters_Int_and_Fire.tau_e
tau_i         = Parameters_Int_and_Fire.tau_i
firing_rate_e = Parameters_Int_and_Fire.firing_rate_e
firing_rate_i = Parameters_Int_and_Fire.firing_rate_i
w_e           = Parameters_Int_and_Fire.w_e
w_i           = Parameters_Int_and_Fire.w_i
delta_t       = Parameters_Int_and_Fire.delta_t

def getPeriodicSpikes(firing_rate,time_step_sim,t_max):
    
    ''' Integer-based method for getting periodic spikes.
    returns desired spike times in step numbers.
    the reason for this method is that a clock with floats will
    lead to skipped spikes due to rounding issues with floats in python.'''
    
    isi_ms = 1000/firing_rate # ISI in ms
    nsteps = int(np.ceil(isi_ms/time_step_sim))
    spikes = np.arange(nsteps,t_max+nsteps,nsteps)
    spiketrue = np.ones([len(spikes)])
    return spikes, spiketrue

# define function for excitatory synapse
def excit_syn(g_e_tt,tau_e): # g_e_tt is the changing variable
    func_g_e=-g_e_tt/tau_e
    return func_g_e
    
# define function for inhibitory synapse
def inhib_syn(g_i_tt,tau_i):
    func_g_i=-g_i_tt/tau_i
    return func_g_i
    
# define function to integrate the membrane voltage equation
arg_for_V='total_input_tt,' + 'E_leak,' + 'tau_mem'
def memb_volt(V_tt,arg_for_V):
    func_V=(E_leak-V_tt+total_input_tt)/tau_mem
    return func_V

# value storage
t_vals = [] # time
V_vals = [] # membrane potential
spike_times = []

# initialize values for voltage, time and conductances (excitatory and inhibitory)
t_vals.append(t_0) 
V_vals.append(V_reset)
V_tt = V_reset
tt = t_0+time_step_sim # initialize tt for the while loop, starts at first time step
g_e_tt_vec=[0]*numb_exc_syn # starting values of all the excitatory conductances (zeros)
g_i_tt_vec=[0]*numb_inh_syn # starting values of all the inhibitory conductances (zeros)
total_exc_input_tt = 0
total_inh_input_tt = 0
buffer_exc_input_tt = []
buffer_inh_input_tt = []
w_exc=[w_e]*numb_exc_syn # define exc. weights (all weights identical)
w_inh=[w_i]*numb_inh_syn # define inh. weights
    
# These counters count the number of excitatory/inhibitory input spikes to the neuron
counter_e = 0
counter_i = 0
number_spikes = 0 # output spike counter

spikes_e, spiketrue_e = getPeriodicSpikes(firing_rate_e,time_step_sim,t_max)
spikes_i, spiketrue_i = getPeriodicSpikes(firing_rate_i,time_step_sim,t_max)

V_1 = Euler()
inputs_e = Input_Regular()
inputs_i = Input_Regular()

###################
# Loop over time
###################

while tt <= t_max:
             
    # create excitatory inputs
    [g_e_tt_vec,counter_e] = inputs_e.inputs(g_e_tt_vec,time_step_sim,spikes_e,spiketrue_e,numb_exc_syn,w_exc,excit_syn,tau_e,tt,counter_e,delta_t)
    total_exc_input_tt = sum([ww*(E_e-V_tt) for ww in g_e_tt_vec]) # get sum of excitatory inputs
    
    # create inhibitory inputs    
    [g_i_tt_vec,counter_i] = inputs_i.inputs(g_i_tt_vec,time_step_sim,spikes_i,spiketrue_i,numb_inh_syn,w_inh,inhib_syn,tau_i,tt,counter_i,delta_t)
    total_inh_input_tt = sum([ww*(E_i-V_tt) for ww in g_i_tt_vec]) # get sum of inhibitory inputs
    
    total_input_tt=total_exc_input_tt+total_inh_input_tt

    V = V_1.euler_integration(memb_volt, arg_for_V, V_tt,tt,time_step_sim,delta_t)
    if V<V_thresh:
        V_tt=V
    else: # if threshold is reached -> set back to V_reset -> spike
        V_tt=V_reset
        V_vals.append(0) # add the spike graphically to the figure
        t_vals.append(tt-time_step_sim/10) # add the spike graphically to the figure
        spike_times.append(tt)
        number_spikes += 1
    
    buffer_exc_input_tt.append(total_exc_input_tt)
    buffer_inh_input_tt.append(total_inh_input_tt)
        
    t_vals.append(tt)
    V_vals.append(V_tt)
    tt=tt+time_step_sim    

##################################
### membrane voltage & conductance
##################################

fig1 = plt.figure()

ax1 = fig1.add_subplot(2,1,1)
plt.plot(t_vals,V_vals)
ax1.set_xlim([0, t_max])
ax1.set_ylabel('Mem. Potential (mV)')
ax1.set_xlabel('Time (ms)')

ax2 = fig1.add_subplot(2,1,2)
ax2.plot(buffer_exc_input_tt,'g', label='Excitation')
ax2.plot(buffer_inh_input_tt,'m',label='Inhibition')
ax2.set_ylabel('Conductance')
ax2.set_xlabel('Time (ms)')
ax2.legend()
plt.tight_layout()
plt.show()
#fig1.savefig('1_2_synaptic_inputs.png')

