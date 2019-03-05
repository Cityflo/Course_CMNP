# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:33:59 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

This code runs the LIF-neuron with a refractory period.
The inputs are excitatory and inhibitory synapses
with Poisson spiking.
Multiple trials are run, and the distributions of ISIs
ad the CV of ISIs over all these trials are shown.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.optimize import curve_fit

from Euler_Method import Euler # import Euler to perform the Euler integration method
from Input_Poisson import Input_Poisson # import class where Poisson input is created

########################################################
# get parameters from parameter file
########################################################

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
E_k           = Parameters_Int_and_Fire.E_k
Delta_g_RP    = Parameters_Int_and_Fire.Delta_g_RP
tau_RP        = Parameters_Int_and_Fire.tau_RP
ntrials       = 100


# define function for excitatory synapse
def excit_syn(g_e_tt,tau_e): # g_e_tt is the changing variable
    func_g_e=-g_e_tt/tau_e
    return func_g_e
    
# define function for inhibitory synapse
def inhib_syn(g_i_tt,tau_i):
    func_g_i=-g_i_tt/tau_i
    return func_g_i
    
# define function to integrate the membrane voltage equation
arg_for_V='total_input_tt,' + 'E_leak,' + 'tau_mem' + 'g_RP_tt' + 'E_k'
def memb_volt(V_tt,arg_for_V):
    func_V=(E_leak-V_tt+total_input_tt+g_RP_tt*(E_k-V_tt))/tau_mem
    return func_V


V_1 = Euler()
inputs_e = Input_Poisson()
inputs_i = Input_Poisson()

CV_all   = []
interspike_intervals_all = []

for trial in range(ntrials):
    
    # value storage
    t_vals = [] # for time
    V_vals = [] # for membrane potential V
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
    w_exc=[w_e]*numb_exc_syn # define exc. weights (all weights identical)
    w_inh=[w_i]*numb_inh_syn # define inh. weights

    # These counters count the number of excitatory/inhibitory input spikes to the neuron
    counter_e=0
    counter_i=0
    number_spikes=0 # counter for the number of output spikes
    g_RP_tt=0
    
    ######################
    # Loop over timepoints
    ######################
    
    while tt <= t_max:
        
        # create excitatory inputs
        [g_e_tt_vec,counter_e]=inputs_e.inputs_pois(g_e_tt_vec,time_step_sim,firing_rate_e,numb_exc_syn,w_exc,excit_syn,tau_e,tt,counter_e,delta_t)
        total_exc_input_tt = sum([ww*(E_e-V_tt) for ww in g_e_tt_vec]) # get sum of excitatory inputs
    
        # create inhibitory inputs
        [g_i_tt_vec,counter_i]=inputs_i.inputs_pois(g_i_tt_vec,time_step_sim,firing_rate_i,numb_inh_syn,w_inh,inhib_syn,tau_i,tt,counter_i,delta_t)
        total_inh_input_tt = sum([ww*(E_i-V_tt) for ww in g_i_tt_vec]) # get sum of inhibitory inputs
    
        total_input_tt=total_exc_input_tt+total_inh_input_tt

        
        V = V_1.euler_integration(memb_volt, arg_for_V, V_tt,tt,time_step_sim,delta_t)
        if V<V_thresh:
            V_tt=V
            g_RP_tt=g_RP_tt-g_RP_tt/tau_RP
        else: # if threshold is reached -> set back to V_reset -> spike
            V_tt = V_reset
            g_RP_tt = g_RP_tt-g_RP_tt/tau_RP+Delta_g_RP
            V_vals.append(0) # add the spike graphically to the figure
            t_vals.append(tt-time_step_sim/10) # add the spike graphically to the figure
            spike_times.append(tt)
            number_spikes = number_spikes+1
        
        if g_RP_tt<0: # g_RP_tt should stay positive
            g_RP_tt = 0
        
        #t_vals.append(tt)
        #V_vals.append(V_tt)
    
        tt=tt+time_step_sim    
        

    ###########################
    # get interspike intervals:
    ###########################
    
    interspike_intervals = []
    interspike_intervals.append(spike_times[0]) # include first interval
    
    for mm in range(1,len(spike_times)):
        interspike_intervals.append(spike_times[mm]-spike_times[mm-1])
        
    interspike_mean = np.mean(interspike_intervals)
    interspike_std = np.std(interspike_intervals)
    coeff_of_variation = interspike_std/interspike_mean
    
    CV_all.append(coeff_of_variation)
    interspike_intervals_all += interspike_intervals
    
#############################################################
# ISI distribution and CV distribution Plots over all trials.
# At least 30 trials are needed to get nice histograms.
############################################################

# fit an exponential function.
def exponential_func(x, a, b):
    return a*np.exp(-b*x)

res = 10
bn = np.arange(0,800+res,res)
bn_plot = np.arange(0,800+10,10)

hst = np.histogram(interspike_intervals_all, bins=bn)
counts = hst[0]
xvals_hist = bn[0:-1] + (bn[1]-bn[0])*0.5
popt, pcov = curve_fit(exponential_func, xvals_hist, counts, p0=[1,0.01])
fit_y = exponential_func(bn_plot, *popt)
print(popt)

fig2 = plt.figure()
st = fig2.suptitle("#Spikes: %(number_spikes)d, Mean ISI: %(interspike_mean)f ms, Mean CV: %(coeff_of_variation)f" % {"number_spikes": number_spikes, "interspike_mean": interspike_mean, "coeff_of_variation": coeff_of_variation})

ax = plt.subplot(1,2,1)

ax.plot(bn_plot,fit_y,lw=3,color='c',label='exp.fit')
ax.hist(interspike_intervals_all,bins=bn,ec ='r',fc = 'None',lw=3,density=False,histtype='step',label='ISIs')
ax.set_xlabel('ISI')
ax.set_ylabel('Frequency')
ax.legend()

ax2 = plt.subplot(1,2,2)
ax2.hist(CV_all,bins=10,lw = 3,histtype='step')
ax2.set_xlabel('CV of ISI')
ax2.set_ylabel('Frequency')
plt.tight_layout()
fig2.subplots_adjust(top=0.90)
plt.show()
fig2.savefig('2_2_RP.png')
