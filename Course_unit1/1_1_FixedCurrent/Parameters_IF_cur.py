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

# simulation time values for the leaky integrate-and-fire neuron:
t_0           = 0 # start time
t_max         = 200 # total simulationtime in ms
time_step_sim = 0.1 # time step for the simulation in ms

# values for the leaky integrate-and-fire neruon:
tau_mem       = 20 # membrane time constant in ms
E_leak        = -60 # reversal potential for the leak in mV
V_reset       = -70 # reset value of the neuron in mV
V_thresh      = -50 # threshold of the neuron in mV
Rm            = 10.0 # resistance in MegaOhm 
Iext          = 2.0 #1.8 # Input current in nanoAmp

# values for the Euler integration
delta_t       = 0.01 # integration step, in ms

# make a file with all parameter values, the date and time
import datetime
now = datetime.datetime.now()

file=open("Parameters.txt","w")
file.write("Parameters used in this simulation:\n")
file.write(" \n") # empty line
file.write("Date: " + str(now) + "\n")
file.write(" \n") # empty line
file.write("tau_mem= %(tau_mem)f\n" % {"tau_mem":tau_mem})
file.write("E_leak= %(E_leak)f\n" % {"E_leak":E_leak})
file.write("V_reset= %(V_reset)f\n" % {"V_reset":V_reset})
file.write("V_thresh= %(V_thresh)f\n" % {"V_thresh":V_thresh})
file.write("Rm= %(Rm)f\n" % {"Rm":Rm})
file.write("Iext= %(Iext)f\n" % {"Iext":Iext})
file.write("t_0= %(t_0)f\n" % {"t_0":t_0})
file.write("t_max= %(t_max)f\n" % {"t_max":t_max})
file.write("time_step_sim= %(time_step_sim)f\n" % {"time_step_sim":time_step_sim})
file.write("delta_t= %(delta_t)f\n" % {"delta_t":delta_t})
file.close()