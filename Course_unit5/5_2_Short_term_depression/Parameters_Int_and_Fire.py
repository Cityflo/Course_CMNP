# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:25:24 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg
"""

#####
# to imrpove it: write values in vectors, call the vectors. for writing: loop over vectors

###############
# Parameter File
###############

# simulation time values for the leaky integrate-and-fire neuron:
t_0=0 # initial time
t_max=5000 # total simulationtime in ms
time_step_sim=1 # time step for the simulation in ms

# values for the leaky integrate-and-fire neruon:
tau_mem=20 # membrane time constant in ms
E_leak=-60 # reversal potential for the leak in mV
V_reset=-70 # reset value of the neuron in mV
V_thresh=-50 # threshold of the neuron in mV

# values for the excitatory conductance g_e differential equation:
E_e=0 # reversal potential for excitatory (depolarizing) inputs in mV
tau_e=3 # postsynaptic potential (PSP), in ms
firing_rate_e=20 #15 # firing rate of the poisson excitatory synapse in Hz (per second!!!!)
w_e=0 # initial strength of each excitatory weight
numb_exc_syn=2 # number of excitatory synapses

# values for the inhibitory conductance g_i differential equation:
E_i=-80 # reversal potential for inhibitory inputs in mV
tau_i=5 # postsynaptic potential (PSP), in ms
firing_rate_i=10 # firing rate of the poisson inhibitory synapse in Hz (per second!!!!)
w_i=0 # initial strength of each inhibitory weights
numb_inh_syn=0 # number of inhibitory synapses

# values for the Euler integration
delta_t=0.01 # integration step, in ms

#### Short-term plasticity  ###
### STD dominated:
#tau_f=50 #ms
#tau_d=750 #ms   
#U=0.45    

### STF dominated:
tau_f=750 #ms    # 50 or 750 ms
tau_d=50 #ms   # 50 or 750 ms
U=0.15    # 0.45 or 0.15 
 
w_fixed = 2.5

import datetime
now = datetime.datetime.now() # gives me date and time

file=open("Parameters.txt","w")
file.write("Parameters used in this simulation:\n")
file.write(" \n") # empty line
file.write("Date: " + str(now) + "\n")
file.write(" \n") # empty line
file.write("tau_mem= %(tau_mem)f\n" % {"tau_mem":tau_mem})
file.write("E_leak= %(E_leak)f\n" % {"E_leak":E_leak})
file.write("V_reset= %(V_reset)f\n" % {"V_reset":V_reset})
file.write("V_thresh= %(V_thresh)f\n" % {"V_thresh":V_thresh})
file.write("E_e= %(E_e)f\n" % {"E_e":E_e})
file.write("tau_e= %(tau_e)f\n" % {"tau_e":tau_e})
file.write("firing_rate_e= %(firing_rate_e)f\n" % {"firing_rate_e":firing_rate_e})
file.write("w_e= %(w_e)f\n" % {"w_e":w_e})
file.write("numb_exc_syn= %(numb_exc_syn)f\n" % {"numb_exc_syn":numb_exc_syn})
file.write("E_i= %(E_i)f\n" % {"E_i":E_i})
file.write("tau_i= %(tau_i)f\n" % {"tau_i":tau_i})
file.write("firing_rate_i= %(firing_rate_i)f\n" % {"firing_rate_i":firing_rate_i})
file.write("w_i= %(w_i)f\n" % {"w_i":w_i})
file.write("numb_inh_syn= %(numb_inh_syn)f\n" % {"numb_inh_syn":numb_inh_syn})
file.write("t_0= %(t_0)f\n" % {"t_0":t_0})
file.write("t_max= %(t_max)f\n" % {"t_max":t_max})
file.write("time_step_sim= %(time_step_sim)f\n" % {"time_step_sim":time_step_sim})
file.write("delta_t= %(delta_t)f\n" % {"delta_t":delta_t})
file.write(" \n") # empty line
file.write("tau_f= %(tau_f)f\n" % {"tau_f":tau_f})
file.write("tau_d= %(tau_d)f\n" % {"tau_d":tau_d})
file.write("U= %(U)f\n" % {"U":U})
file.write("w_fixed= %(w_fixed)f\n" % {"w_fixed":w_fixed})
file.close()