import numpy as np
from matplotlib.pyplot import figure, plot, show, imshow
import matplotlib.pyplot as plt
import os

os.chdir('/Users/danhobley/Documents/SIESD/delta_expt/overhead images/all_data')

my_tslices = np.arange(250, 5650, 20, dtype=int)

# get the first one
first_one = np.loadtxt('800_not_in_channel_'+str(my_tslices[0])+'.csv', delimiter=',')

# make a stencil
stencil = np.ones_like(first_one, dtype=bool)
channel_prop_record = []
not_subaerial_prop_record = []


# kill the sides
stencil[:100,:] = False
stencil[:,:100] = False

# slice a diagonal stripe
for i in xrange(400):
    stencil[i,:(400-i)] = False
    #stencil[]
for i in xrange(600):
    stencil[i,(600-i):] = False
stencil[600:,:] = False

stencil_size = stencil.sum()

for mytslice in my_tslices:
    print mytslice
    not_channel = np.loadtxt('800_not_in_channel_'+str(mytslice)+'.csv', delimiter=',')
    subaerial = np.loadtxt('800_subaerial_'+str(mytslice)+'.csv', delimiter=',')
    channel_cells = np.logical_and(np.logical_not(not_channel), stencil)
    not_subaerial_cells = np.logical_and(np.logical_not(subaerial), stencil)
    channel_prop = float(channel_cells.sum())/stencil_size
    not_subaerial_prop = float(not_subaerial_cells.sum())/stencil_size
    channel_prop_record.append(float(channel_prop))
    not_subaerial_prop_record.append(float(not_subaerial_prop))

channel_prop_record = np.array(channel_prop_record)
not_subaerial_prop_record = np.array(not_subaerial_prop_record)

channel_prop_median = np.median(channel_prop_record)
not_subaerial_prop_median = np.median(not_subaerial_prop_record)

os.chdir('/Users/danhobley/stochastic-delta')
