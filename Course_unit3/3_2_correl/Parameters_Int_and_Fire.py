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
t_max=15000 # total simulationtime in ms
time_step_sim=1 # time step for the simulation in ms

# values for the excitatory inputs:
firing_rate_e=10 # firing rate of the poisson excitatory synapse in Hz (per second!!!!)
numb_exc_syn=20 # number of excitatory synapses

# values for the correlated spike trains
c1 = 0.2
c2 = 0.1
tau_c = 20 # ms

import datetime
now = datetime.datetime.now() # gives me date and time

file=open("Parameters.txt","w")
file.write("Parameters used in this simulation:\n")
file.write(" \n") # empty line
file.write("Date: " + str(now) + "\n")
file.write(" \n") # empty line
file.write("firing_rate_e= %(firing_rate_e)f\n" % {"firing_rate_e":firing_rate_e})
file.write("numb_exc_syn= %(numb_exc_syn)f\n" % {"numb_exc_syn":numb_exc_syn})
file.write("t_0= %(t_0)f\n" % {"t_0":t_0})
file.write("t_max= %(t_max)f\n" % {"t_max":t_max})
file.write("time_step_sim= %(time_step_sim)f\n" % {"time_step_sim":time_step_sim})
file.write(" \n") # empty line
file.write("c1= %(c1)f\n" % {"c1":c1})
file.write("c2= %(c2)f\n" % {"c2":c2})
file.write("tau_c= %(tau_c)f\n" % {"tau_c":tau_c})
file.close()