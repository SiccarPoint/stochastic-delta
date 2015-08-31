# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:50:49 2015

Simple script to run the delta.

@author: danhobley
"""

import numpy as np
import delta as d
from pylab import show, plot, figure
import matplotlib.pyplot as plt

#history = np.array([1.,.9,.8,.7,.6])
#history = np.array([1.,1.])
#history = np.array([1.,.9,.8])
#history = np.array([2.,0.5])
#history = np.array([0.5,2.])
history = np.random.rand(100)*2.+0.2
#history = np.array([1.,1.,0.5,2.,0.1,0.6])
#history = np.array([0.90275123,  0.77150464,  0.81842169])
#history = np.array([0.77150464,  0.81842169])
#history = np.ones(100, dtype=float)*8.
#history[4:] = 1.  #agressive SL drop on existing structure
mydelta = d.delta_shoreline(history)
mydelta.execute(500.)

figure(3)
for i in xrange(history.size):
    plot((mydelta.toe_x[i], mydelta.all_x[i], mydelta.head_x[i]),
         (mydelta.toe_z[i], mydelta.all_z[i], mydelta.head_z[i]),'-')

figure(4)
for i in xrange(history.size):
    if mydelta.strat_foresets_x[0,i] != -1.:
        plot(mydelta.strat_foresets_x[:,i], mydelta.strat_foresets_z[:,i], 'k-')
    if mydelta.strat_topsets_x[0,i] != -1.:
        plot(mydelta.strat_topsets_x[:,i], mydelta.strat_topsets_z[:,i], 'k-')

figure(5)
for i in xrange(history.size):
    plot((mydelta.topo_toe_x[i], mydelta.all_x[i], mydelta.topo_head_x[i]),
         (mydelta.topo_toe_z[i], mydelta.all_z[i], mydelta.topo_head_z[i]),'-')

mydelta.full_completeness()
