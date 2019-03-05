# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:33:59 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

This code runs a simple LIF-neuron model with fixed current input.
"""


# import Euler to perform the Euler integration method
from Euler_Method import Euler
import Parameters_IF_cur

########################################################

########################################################

tau_mem       = Parameters_IF_cur.tau_mem
E_leak        = Parameters_IF_cur.E_leak
V_reset       = Parameters_IF_cur.V_reset
V_thresh      = Parameters_IF_cur.V_thresh
t_0           = Parameters_IF_cur.t_0
t_max         = Parameters_IF_cur.t_max
time_step_sim = Parameters_IF_cur.time_step_sim
Rm            = Parameters_IF_cur.Rm
Iext          = Parameters_IF_cur.Iext 
delta_t       = Parameters_IF_cur.delta_t

# define function to integrate the membrane voltage equation
arg_for_V = 'E_leak,' + 'tau_mem' + 'Rm' + 'Iext'
def memb_volt(V_tt,arg_for_V):
    func_V=(E_leak-V_tt+(Rm*Iext))/tau_mem
    return func_V

# store the values
x_vals=[]
y_vals=[]
x_vals.append(t_0) # add first element
y_vals.append(V_reset) # add first element

# initialize values for voltage, time and excitatory conductance
V_tt=V_reset # start value of the voltage
tt=t_0+time_step_sim # initialize tt for the while loop, starts at first time step
number_spikes=0 # counter for the number of spikes

###################
# Loop over time
###################

# call function from Euler_Method to integrate membrane voltage equation
V_1=Euler()

# a while loop also allows non-integer step sizes.
while tt <= t_max:
    
    if V_tt==0: # if spike -> set back to V reset
        V_tt=V_reset
        x_vals.append(tt)
        y_vals.append(V_tt)
    else:       
        V=V_1.euler_integration(memb_volt, arg_for_V, V_tt,tt,time_step_sim,delta_t)
        
        # if threshold is reached -> set back to V_reset -> spike
        if V<V_thresh:
            V_tt=V
        else:
            V_tt=0 # a spike is shown with peak value V = 0 mV.
            print('We have a spike')
            number_spikes=number_spikes+1
        
        x_vals.append(tt) # storage of time values
        y_vals.append(V_tt) # storage of V_tt values
        
        tt=tt+time_step_sim # proceed to next time step

print('Number of spikes:')
print(number_spikes)

###################
# Plot membrane potential
###################

import matplotlib.pyplot as plt

fig1 = plt.figure()
plt.plot(x_vals,y_vals)
plt.ylabel('Mem. Potential (mV)')
plt.xlabel('Time (ms)')
plt.show()
fig1.savefig('1_1_current_input.png')

