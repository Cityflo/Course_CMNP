# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:33:59 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

This simulation contains a LIF neuron with STDP and
short-term plasticity (STP) in its excitatory input synapses.
Here, the STP includes only short-term facilitation (STF).
Inhibitory inputs are omitted for simplicity.

The STP model is based on:
Tsodyks, M. et al. (1998), Neural Networks with Dynamic Synapses. Neural Computation 10(4): 821-835

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

from Euler_Method import Euler 
from Input_synapse import Synapse

from Poisson_Spike_Trains import Poisson_Trains

from Plasticity_STP import STPlasticity

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
# STP parameters
tau_f         = Parameters_Int_and_Fire.tau_f
tau_d         = Parameters_Int_and_Fire.tau_d
U             = Parameters_Int_and_Fire.U
w_fixed       = Parameters_Int_and_Fire.w_fixed

# define function for excitatory synapse
def excit_cond(g_e_tt,tau_e): # g_e_tt is the changing variable
    func_g_e=-g_e_tt/tau_e
    return func_g_e

# define function for inhibitory synapse
def inhib_cond(g_i_tt,tau_i): # g_e_tt is the changing variable
    func_g_i=-g_i_tt/tau_i
    return func_g_i

# define function to integrate the membrane voltage equation
arg_for_V='total_input_tt,' + 'E_leak,' + 'tau_mem' + 'E_k'
def memb_volt(V_tt,arg_for_V):
    func_V=(E_leak-V_tt+total_input_tt)/tau_mem
    return func_V

# value storage
t_vals = [] # for time
V_vals = [] # for membrane potential V
spike_times = []
FR_vec=[0]

# initialize values for voltage, time and conductances (excitatory and inhibitory)
t_vals.append(t_0) 
V_vals.append(V_reset)
V_tt = V_reset
tt = t_0+time_step_sim # initialize tt for the while loop, starts at first time step
g_e_tt_vec=[0]*numb_exc_syn # starting values of all the excitatory conductances (zeros)
g_i_tt_vec=[0]*numb_inh_syn # starting values of all the inhibitory conductances (zeros)
total_exc_input_tt = 0
total_inh_input_tt = 0
w_e_vec_tt=[w_e]*numb_exc_syn # define exc. weights here (all the weights are the same)
w_i_vec_tt=[w_i]*numb_inh_syn 

# buffer for spike times (for STDP)
last_post_spike=0
idx_of_synapse_with_spike=[]

# buffer for weight changes over time (STDP in excitatory synapses)
w_e_storage=np.zeros((int(round((t_max-t_0)/time_step_sim))+1, numb_exc_syn))
w_e_storage[0,:]=w_e_vec_tt
counter_storage=1;

# These counters count the number of excitatory/inhibitory inputs to the neuron
counter_e=0
spike_time_e=[0]*numb_exc_syn
counter_i=0
spike_time_i=[0]*numb_inh_syn
number_spikes = 0 

# For STP:
spike_or_not=0 # used for STP. did a pre-spike happen?
u_vec=[0]*numb_exc_syn
x_vec=[1]*numb_exc_syn # x remains fixed in this simulation. (no STD)
u_storage=np.zeros((int(round((t_max-t_0)/time_step_sim))+1, numb_exc_syn))
u_storage[0,:]=u_vec

###########################
# create input spike trains
###########################

# firing rates : 
r1 = firing_rate_e
r2 = 0#firing_rate_e # silence other synapse

spikes_e = Poisson_Trains()
[list_of_all_spike_trains1,list_of_all_spike_trains2] = spikes_e.get_list_of_trains(r1,r2)
spike_trains_complete_e = list_of_all_spike_trains1 + list_of_all_spike_trains2

###################
# Start of the loop
###################

synapse_e = Synapse()

exc_STP = STPlasticity()

# use while to make sure that also step sizes < 1 are possible (that's a problem in python for loops)
while tt <= t_max:

    # call function for excitatory inputs
    [g_e_tt_vec]=synapse_e.inputs_calc(g_e_tt_vec,time_step_sim,numb_exc_syn,w_e_vec_tt,excit_cond,tau_e,tt,delta_t,spike_trains_complete_e)
    total_exc_input_tt=sum([ww*(E_e-V_tt) for ww in g_e_tt_vec])
    
    # Include STP in the excitatory synapses.
    [w_e_vec_tt,u_vec,x_vec] = exc_STP.STP(spike_trains_complete_e,tt,time_step_sim,tau_f,tau_d,U,w_fixed,w_e_vec_tt,u_vec,x_vec)
        
    # update the postsynaptic neuron.
    total_input_tt=total_exc_input_tt
    V_1 = Euler()
    V = V_1.euler_integration(memb_volt, arg_for_V, V_tt,tt,time_step_sim,delta_t)
    
    if V<V_thresh:
        V_tt=V
    else: # if threshold is reached -> spike & reset mem V
        V_tt=V_reset
        V_vals.append(0) 
        t_vals.append(tt-time_step_sim/10) 
        spike_times.append(tt)
        number_spikes=number_spikes+1

    t_vals.append(tt)
    V_vals.append(V_tt)
    tt=tt+time_step_sim
    
    u_storage[counter_storage,:] = u_vec
    spike_or_not = 0 # set back to "no spike"
    
    # store firing rate every second
    if tt%1000==0:
        R_count=number_spikes/tt*1000
        FR_vec.append(R_count)
    
    w_e_storage[counter_storage,:]=w_e_vec_tt # store the weights of the excitatory synapses
    counter_storage=counter_storage+1

######################################################
# Plot mem. potential and weight evolution
######################################################

ex_syn  = 0  # choose a synapse.

fig1=plt.figure()
gs = gridspec.GridSpec(3, 1)
# show dynamics for the single synapse weight
ax0 = fig1.add_subplot(gs[0,0])
ax0.plot(range(int(round((t_max-t_0)/time_step_sim))+1),u_storage[:,ex_syn],lw=2,color='k')
ax0.set_xlabel('Time (ms)')
ax0.set_ylabel('u')

ax1 = fig1.add_subplot(gs[1,0])
ax1.plot(range(int(round((t_max-t_0)/time_step_sim))+1),w_e_storage[:,ex_syn],lw=2)
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Syn. Weight')

# postsynaptic firing rate
ax2 = fig1.add_subplot(gs[2,0])
ax2.plot(FR_vec,color='r',lw=2)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Firing Rate (Hz)')

plt.tight_layout()
plt.show()
fig1.savefig('5_1_STF_'+ str(tau_f) + 'ms.png')
