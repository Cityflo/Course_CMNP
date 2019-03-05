# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:32:43 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

"""
#!/usr/bin/python

################
# Description:
# Class that performs the Euler integration. Call the function euler_integration 
# and pass it with the function f which is dy/dt=f(y,t), the initial value y_0 and 
# t_0 and the total time for which the integration should be done (time_step_sim).
# The output of the function will be the value after the time time_step_sim
# and the total time t_1 at this value.
################

class Euler: 

    def euler_integration(self, f, arg_for_f, y_0, t_0, time_step_sim, delta_t):
        
        if delta_t>=time_step_sim:
            a=int(input("ATTENTION: The time step of the simulation is smaller than the time step of the integration!"))
            
        n=int(round(time_step_sim/delta_t)) # calculate the number of steps
        for step in range(1,n+1):
            m=f(y_0, arg_for_f)
            y_1=y_0+delta_t*m
            t_1=t_0+delta_t
            t_0=t_1
            y_0=y_1

        return y_0