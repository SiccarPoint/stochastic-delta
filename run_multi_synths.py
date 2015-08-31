from delta_obj2 import delta
import numpy as np
from matplotlib.pyplot import colorbar, figure, show, plot, imshow
import matplotlib.pyplot as plt

num_versions = 30

mydelta = delta()
ins = mydelta.read_input_file('synth_inputs.txt')
nt = int(ins['nt'])
SL_trajectory = (np.random.rand(nt)*2.+0.2)
SL_trajectory += np.arange(nt)*0.01
SL_trajectory2 = np.sin(np.arange(nt)/2.)+1.2
SL_trajectory2 += np.arange(nt)*0.01
SL_trajectory2 += np.sin(np.arange(nt)/4.)*0.5
SL_trajectory3 = np.sin(np.arange(nt)/2.3)
SL_trajectory3 *= np.sin(np.arange(nt)/2.1)
SL_trajectory3 *= np.sin(np.arange(nt)/4.1)*0.5
SL_trajectory3 *= 2.
SL_trajectory3 += np.arange(nt)*0.01 + 0.7
completenesses_walking = []
completenesses_whole = []
completenesses_restricted = []
completenesses_noE = []

for i in xrange(num_versions):
    tscales, completenesses_walking = mydelta.execute('test_inputs.txt', SL_trajectory3, completeness_records=completenesses_walking, graphs=False, walking_erosion_depo=True)
    tscales, completenesses_restricted = mydelta.execute('test_inputs.txt', SL_trajectory3, completeness_records=completenesses_restricted, graphs=False, restricted_channel_mass_conserved=True)
tscales, completenesses_whole = mydelta.execute('test_inputs.txt', SL_trajectory3, completeness_records=completenesses_whole, graphs=False)
tscales, completenesses_noE = mydelta.execute('test_inputs.txt', SL_trajectory3, graphs=False, completeness_records=completenesses_noE, never_erosion=True)
mean_comp_walking = np.mean(completenesses_walking, axis=0)
mean_comp_whole = np.mean(completenesses_whole, axis=0)
mean_comp_restricted = np.mean(completenesses_restricted, axis=0)
figure(1)
plt.gca().set_xscale('log')
plt.xlabel('multiple of smallest tstep')
plt.ylabel('completeness')
for i in completenesses_walking:
    plot(tscales, i)
figure(2)
plt.gca().set_xscale('log')
plt.xlabel('multiple of smallest tstep')
plt.ylabel('completeness')
for i in completenesses_restricted:
    plot(tscales, i)
figure(3)
plt.gca().set_xscale('log')
plt.xlabel('multiple of smallest tstep')
plt.ylabel('completeness')
for i in completenesses_whole:
    plot(tscales, i)
figure(4)
plt.gca().set_xscale('log')
plt.xlabel('multiple of smallest tstep')
plt.ylabel('completeness')
plot(tscales, completenesses_noE[0])
figure(5)
plt.gca().set_xscale('log')
plt.xlabel('multiple of smallest tstep')
plt.ylabel('completeness')
plot(tscales, mean_comp_walking)
plot(tscales, mean_comp_whole)
plot(tscales, mean_comp_restricted)
plot(tscales, completenesses_noE[0])

figure(6)
plt.gca().set_xscale('log')
plt.xlabel('multiple of smallest tstep')
plt.ylabel('completeness')
for i in completenesses_walking:
    plot(tscales, i, color=(0.35,0.35,1.), ls='-.', lw=1.)
for i in completenesses_restricted:
    plot(tscales, i, color=(0.45,1.,0.45), ls='--', lw=1.)
plot(tscales, mean_comp_walking, color=(0.,0.,1.), lw=3.)
plot(tscales, mean_comp_whole, color=(1.,0.,0.), lw=3.)
plot(tscales, mean_comp_restricted, color=(0.,1.,0.), lw=3.)
plot(tscales, completenesses_noE[0], color='k', lw=3.)
