# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 12:33:59 2017
Tested in Python 3.6.5, Anaconda Inc.

@author: Miehl
Adapted by: Florence Kleberg

This code creates separate groups of correlated
spike trains. The LIF model is not included here.

Correlation within a group is identical, 
but can be different beween groups. 
Correlated_Trains() is used for instantaneous correlations.
CorrelatedJitter_Trains() is used for correlations with exponential 
jitter, leading to broader cross-correlograms.

The cross-correlation a pair of spike trains, from either the same or a different group, 
is shown in the figure, where the zero-lag is indicated
by a vertical line.
"""

import numpy as np
import matplotlib.pyplot as plt

### create correlated spike trains for excitatory inputs.
from Correlated_Spike_Trains import Correlated_Trains
from CorrelatedJitter_Spike_Trains import CorrelatedJitter_Trains

import Parameters_Int_and_Fire

t_max         = Parameters_Int_and_Fire.t_max
firing_rate_e = Parameters_Int_and_Fire.firing_rate_e

# correlation in the two groups
c1            = Parameters_Int_and_Fire.c1     
c2            = Parameters_Int_and_Fire.c2
tau_c         = Parameters_Int_and_Fire.tau_c

#### get correlated spike trains:

### instantaneous correlations
#spikes_e_corr = Correlated_Trains()
#[list_of_all_spike_trains1,list_of_all_spike_trains2] = spikes_e_corr.get_list_of_trains(c1,c2,firing_rate_e)

### jittered correlations
spikes_e_corr = CorrelatedJitter_Trains()
[list_of_all_spike_trains1,list_of_all_spike_trains2] = spikes_e_corr.get_list_of_trains(c1,c2,firing_rate_e,tau_c)

###############################################

spike_trains_complete_e = list_of_all_spike_trains1 + list_of_all_spike_trains2

###############################################
# bin spikes (in order to get correlations)
# Note that correlations become larger
# as binsize increases.
###############################################

binning_bin = 5 # in ms
# decide on which spike trains are cross-correlated:
# same group or different group? See parameter file for group size.
train1 = 0
train2 = 1 

def BinSpikesFast(spiketrain,binsize,run_time):
    # loopless version for binning spikes. Requires numpy
    binlabels = np.arange(0,int(run_time)+binsize,binsize)
    binvalues = np.zeros([len(binlabels)],float)
        
    mat = np.tile(binlabels, (len(spiketrain),1))    
    mat_err = np.absolute(mat - np.transpose(np.tile(spiketrain,(len(binlabels),1))))
    
    # find index of bin with smallest error. If two neighouring bins have equal error,
    # then the numpy.argmin() function always chooses the lowest bin.
    bin_ind = np.argmin(mat_err,axis=1)

    # some spikes end up in the same bin.
    # to avoid loop argument, 
    # use numpy.histogram to find total spikes that go in each bin.
    # the histogram classifies just integers (the indices)
    hst = np.histogram(bin_ind, bins = binlabels)
    return hst[0] # how many spikes in each bin.

binned_spikes1 = BinSpikesFast(np.array(spike_trains_complete_e[train1]),binning_bin,t_max)
binned_spikes2 = BinSpikesFast(np.array(spike_trains_complete_e[train2]),binning_bin,t_max)

##########################################################
# Plot the cross-correlations between the two spike trains
##########################################################

CORR = np.correlate(binned_spikes1,binned_spikes2,'full')

plot_range = 100 # ms
central_idx = int(round((t_max/binning_bin)))-1
plot_xmin = np.floor(central_idx - plot_range/binning_bin)
plot_xmax = np.floor(central_idx + plot_range/binning_bin)

corr_plot = CORR[int(plot_xmin):int(plot_xmax+1)]
corr_plot = corr_plot / (firing_rate_e**2)# normalise by the firing rates

fig1 = plt.figure()
plt.plot(corr_plot,lw=2)
plt.xticks([0,plot_range/binning_bin,2*plot_range/binning_bin],[-plot_range,0,plot_range])
zero_lag = plot_range/binning_bin # CORR[central_idx]/ (firing_rate_e**2)
plt.plot([zero_lag,zero_lag],[min(corr_plot),max(corr_plot)],'--',lw=2,color='m') # zero-lag correlation
plt.xlabel('Cross-correlation lag (ms)')
plt.ylabel('Correlation')
plt.tight_layout()
plt.show()
fig1.savefig('3_2_cross_correlation.png')

############################################################################################################

