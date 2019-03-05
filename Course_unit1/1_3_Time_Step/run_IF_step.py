# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:33:59 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

This code is for testing the accuracy of the Euler integration
and choosing the correct time step.
For this, the exponential synapse model from 1.2 is used.  
"""

import numpy as np
import matplotlib.pyplot as plt
from Euler_Method_Check import Euler
import Parameters_IF_step

tau_e         = Parameters_IF_step.tau_e
g_e_0         = Parameters_IF_step.g_e_0
t_0           = Parameters_IF_step.t_0
t_max         = Parameters_IF_step.t_max
time_step_sim = Parameters_IF_step.time_step_sim
delta_t       = Parameters_IF_step.delta_t
delta_t_bad   = Parameters_IF_step.delta_t_bad

# store the values
x_vals=[]
y_vals=[]

# add initial values
x_vals.append(t_0) 
y_vals.append(g_e_0)

# calculate the results of the exponential function e^(-t/tau_e)
res = 0.0001 # resolution of analytical solution (choose very small!)
analytical_solution_x= np.arange(t_0,t_max,res)
analytical_solution_y= np.exp(-analytical_solution_x/tau_e)

# define function for exponential excitatory synapse
def excit_syn(g_e_tt,tau_e):
    func_g_e=-g_e_tt/tau_e
    return func_g_e
    
#####################
# Loop over timesteps
#####################
y_1=Euler()

tt = t_0
while tt <= t_max:
             
    [g_e,t_1]=y_1.euler_integration(excit_syn, tau_e, g_e_0, tt, time_step_sim, delta_t)

    x_vals.append(t_1)
    y_vals.append(g_e)
    g_e_0 = g_e # use g_e as new input value for the integration
    
    tt = tt+time_step_sim    
    
#####################    
### Show an example of a bad timestep, and reset g_e_0:
#####################

x_vals_bad=[]
y_vals_bad=[]
tt = t_0
g_e_0  = Parameters_IF_step.g_e_0

while tt <= t_max:
             
    [g_e_bad,t_1]=y_1.euler_integration(excit_syn, tau_e, g_e_0, tt, time_step_sim, delta_t_bad)
    
    x_vals_bad.append(t_1)
    y_vals_bad.append(g_e_bad)
    g_e_0 = g_e_bad # use g_e as new input value for the integration
    
    tt = tt+time_step_sim    

#####################
### plot results
#####################

fig1 = plt.figure()
plt.plot(x_vals,y_vals,label='Sim, dt = ' + str(delta_t) + ' ms')
plt.plot(x_vals_bad,y_vals_bad,label='Sim, dt = ' + str(delta_t_bad) + ' ms',color=[0.7,0.7,0.7])
plt.plot(analytical_solution_x,analytical_solution_y,'--',label='Ana. solution')
plt.ylabel('g_e',)
plt.xlabel('time t in ms')
plt.xlim(0,30)
plt.legend()
plt.title('Simulation time step: ' + str(time_step_sim) + ' ms')
plt.tight_layout()
plt.show()
fig1.savefig('1_3_integration_time_step.png')


