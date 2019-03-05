# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:33:59 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

This code runs a LIF-neuron model with Spike-Rate Adaptation (SRA). 
There are no excitatory and inhibitory synapses here.
The input is fixed current input, as in unit 1.1, this makes 
the effect of SRA easier to see. 
The gradual decrease in frequency is shown in the last subfigure.
"""

from Euler_Method import Euler
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import Parameters_IF

tau_mem       = Parameters_IF.tau_mem
E_leak        = Parameters_IF.E_leak
V_reset       = Parameters_IF.V_reset
V_thresh      = Parameters_IF.V_thresh
t_0           = Parameters_IF.t_0
t_max         = Parameters_IF.t_max
time_step_sim = Parameters_IF.time_step_sim
delta_t       = Parameters_IF.delta_t
E_k           = Parameters_IF.E_k
Delta_g_SRA   = Parameters_IF.Delta_g_SRA
tau_SRA       = Parameters_IF.tau_SRA
Rm            = Parameters_IF.Rm
Iext          = Parameters_IF.Iext 
cur_start     = Parameters_IF.cur_start
cur_stop      = Parameters_IF.cur_stop


# define function to integrate the membrane voltage equation
arg_for_V='E_leak,' + 'tau_mem' + 'g_SRA_tt' + 'E_k'
def memb_volt(V_tt,arg_for_V):
    func_V=(E_leak-V_tt+g_SRA_tt*(E_k-V_tt)+(Rm*Iext))/tau_mem
    return func_V

# value storage
t_vals = [] # for time
t_vals_fr = [] # for timing of firing rate recording
V_vals = [] # for membrane potential V
spike_times = []
FR = []
g_SRA_vals = []

# initialize values for voltage, time and conductances (excitatory and inhibitory)
t_vals.append(t_0) 
V_vals.append(V_reset)
g_SRA_vals.append(0)
V_tt = V_reset
tt = t_0+time_step_sim # initialize tt for the while loop, starts at first time step
number_spikes=0 # counter for the number of output spikes
g_SRA_tt=0
    
######################
# Loop over timepoints
######################

while tt <= t_max:

    if tt < cur_start or tt > cur_stop:
        Iext=0
    else:
        Iext = Parameters_IF.Iext   

    V_1 = Euler()
    V = V_1.euler_integration(memb_volt, arg_for_V, V_tt,tt,time_step_sim,delta_t)
    
    if V<V_thresh:
        V_tt = V
        g_SRA_tt=g_SRA_tt-g_SRA_tt/tau_SRA
    else: # if threshold is reached -> set back to V_reset -> spike
        V_tt=V_reset
        # include SRA:
        g_SRA_tt=g_SRA_tt-g_SRA_tt/tau_SRA+Delta_g_SRA
        # extra values for plotting spikes: & consistency with g_SRA_vals
        V_vals.append(0)
        t_vals.append(tt-time_step_sim/10) 
        g_SRA_vals.append(g_SRA_tt)
        spike_times.append(tt)
        number_spikes += 1
    
    t_vals.append(tt)
    V_vals.append(V_tt)
    g_SRA_vals.append(g_SRA_tt)
    tt += time_step_sim    
    
fig1 = plt.figure()
gs = gridspec.GridSpec(2, 2)

ax1 = fig1.add_subplot(gs[0,0])
ax1.plot(t_vals,V_vals)
ax1.plot([cur_start,cur_start],[-70,0],'--',color='k')
ax1.set_xlabel('time (ms)')
ax1.set_ylabel('Mem. potential (mV)')

ax2 = fig1.add_subplot(gs[1,0])
ax2.plot(t_vals,g_SRA_vals)
ax2.set_xlabel('time (ms)')
ax2.set_ylabel('g_SRA')

ax3 = fig1.add_subplot(gs[:,1])
ISIs = np.array(spike_times) - np.array([cur_start] + spike_times[0:-1])
ax3.plot(np.arange(1,len(spike_times)+1,1),1000/ISIs)
ax3.set_xlabel('Spike number')
ax3.set_ylabel('Spike frequency (Hz)')
ax3.set_ylim(20,50)
plt.tight_layout()
plt.show()
if Delta_g_SRA > 0:
    fig1.savefig('2_1_SRA_V.png')
else:
    fig1.savefig('2_1_SRA_V_noSRA.png')
