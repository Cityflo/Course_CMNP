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
t_0=0 # initial time
t_max=500 # total simulationtime in ms
time_step_sim=1 # time step for the simulation in ms
cur_start = 50 # start input current
cur_stop  = 350 # stop input current

# values for the leaky integrate-and-fire neruon:
tau_mem=20 # membrane time constant in ms
E_leak=-60 # reversal potential for the leak in mV
V_reset=-70 # reset value of the neuron in mV
V_thresh=-50 # threshold of the neruon in mV

# values for the spike-rate adoption
E_k=-70 # in mV
Delta_g_SRA= 0.06
tau_SRA=100 #ms
# values for the input current
Rm = 10.0 # resistance in MegaOhm 
Iext = 1.45 # Input current in nanoAmp

# values for the Euler integration
delta_t=0.01 # integration step, in ms

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
file.write("t_0= %(t_0)f\n" % {"t_0":t_0})
file.write("t_max= %(t_max)f\n" % {"t_max":t_max})
file.write("time_step_sim= %(time_step_sim)f\n" % {"time_step_sim":time_step_sim})
file.write("cur_start= %(cur_start)f\n" % {"cur_start":cur_start})
file.write("cur_stop= %(cur_stop)f\n" % {"cur_stop":cur_stop})
file.write("E_k= %(E_k)f\n" % {"E_k":E_k})
file.write("Delta_g_SRA= %(Delta_g_SRA)f\n" % {"Delta_g_SRA":Delta_g_SRA})
file.write("tau_SRA= %(tau_SRA)f\n" % {"tau_SRA":tau_SRA})
file.write("delta_t= %(delta_t)f\n" % {"delta_t":delta_t})
file.write("Rm= %(Rm)f\n" % {"Rm":Rm})
file.write("Iext= %(Iext)f\n" % {"Iext":Iext})
file.close()