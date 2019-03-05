# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:25:24 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

"""

###############
# Parameter File
###############

# values for the leaky integrate-and-fire neuron:
tau_e         = 3 # postsynaptic potential (PSP), in ms
g_e_0         = 1 # value of excitatory conductance at time t=0
t_0           = 0 # initial time
t_max         = 100 # total simulationtime in ms
time_step_sim = 1 # time step for the simulation in ms
delta_t       = 0.01 # Euler integration step, in ms
delta_t_bad   = 1 # Example of badly chosen Euler integration step, in ms

import datetime

now = datetime.datetime.now() # gives me date and time

file=open("Parameters.txt","w")
file.write("Parameters used in this simulation:\n")
file.write(" \n") # empty line
file.write("Date: " + str(now) + "\n")
file.write("tau_e= %(tau_e)f\n" % {"tau_e":tau_e})
file.write("g_e_0= %(g_e_0)f\n" % {"g_e_0":g_e_0})
file.write("t_0= %(t_0)f\n" % {"t_0":t_0})
file.write("t_max= %(t_max)f\n" % {"t_max":t_max})
file.write("time_step_sim= %(time_step_sim)f\n" % {"time_step_sim":time_step_sim})
file.write("delta_t= %(delta_t)f\n" % {"delta_t":delta_t})
file.write("delta_t_bad= %(delta_t_bad)f\n" % {"delta_t_bad":delta_t_bad})
file.close()