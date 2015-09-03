from delta_obj2 import delta
import numpy as np
from matplotlib.pyplot import colorbar, figure, show, plot, imshow, legend, title
import matplotlib.pyplot as plt

num_versions = 10

mydelta = delta()
ins = mydelta.read_input_file('sensitivity_inputs.txt')
nt = int(ins['nt'])
SL_trajectory3 = np.sin(np.arange(nt)/2.3)
SL_trajectory3 *= np.sin(np.arange(nt)/2.1)
SL_trajectory3 *= np.sin(np.arange(nt)/4.1)*0.5
SL_trajectory3 *= 2.
# SL_trajectory3 += np.arange(nt)*0.01 + 0.7
SL_trajectory3 += 2.
SL_miller = (np.loadtxt('Miller_SL_1000y.txt')[:nt])[::-1]+150.
# ^note the reversal of order!!!
completenesses_walking = []
completenesses_whole = []
completenesses_restricted = []
completenesses_noE = []
completenesses_compensated = []

# for i in xrange(num_versions):
#     tscales, completenesses_walking = mydelta.execute('synth_inputs.txt', SL_trajectory3, completeness_records=completenesses_walking, graphs=False, walking_erosion_depo=True)
#     tscales, completenesses_restricted = mydelta.execute('synth_inputs.txt', SL_trajectory3, completeness_records=completenesses_restricted, graphs=False, restricted_channel_mass_conserved=True)
#     tscales, completenesses_compensated = mydelta.execute('synth_inputs.txt', SL_trajectory3, completeness_records=completenesses_compensated, graphs=False, compensation=True)
tscales, completenesses_whole = mydelta.execute('sensitivity_inputs.txt', SL_trajectory=SL_miller, completeness_records=completenesses_whole, graphs=True)
# tscales, completenesses_noE = mydelta.execute('synth_inputs.txt', SL_trajectory3, graphs=False, completeness_records=completenesses_noE, never_erosion=True)
# mean_comp_walking = np.mean(completenesses_walking, axis=0)
mean_comp_whole = np.mean(completenesses_whole, axis=0)
# mean_comp_restricted = np.mean(completenesses_restricted, axis=0)
# mean_comp_compensated = np.mean(completenesses_compensated, axis=0)
gray = (0.7,0.7,0.7)
# figure(1)
# plt.gca().set_xscale('log')
# plt.xlabel('multiple of smallest tstep')
# plt.ylabel('completeness')
# for i in completenesses_walking:
#     plot(tscales, i, color=gray)
# plot(tscales, mean_comp_walking, 'k', lw=3.)
# title('walking')
# figure(2)
# plt.gca().set_xscale('log')
# plt.xlabel('multiple of smallest tstep')
# plt.ylabel('completeness')
# for i in completenesses_restricted:
#     plot(tscales, i, color=gray)
# plot(tscales, mean_comp_restricted, 'k', lw=3.)
# title('restricted')
# figure(3)
# plt.gca().set_xscale('log')
# plt.xlabel('multiple of smallest tstep')
# plt.ylabel('completeness')
# for i in completenesses_whole:
#     plot(tscales, i, color=gray)
# plot(tscales, mean_comp_whole, 'k', lw=3.)
# title('whole')
# figure(4)
# plt.gca().set_xscale('log')
# plt.xlabel('multiple of smallest tstep')
# plt.ylabel('completeness')
# plot(tscales, completenesses_noE[0], 'k', lw=3.)
# title('no erosion')
# figure(5)
# plt.gca().set_xscale('log')
# plt.xlabel('multiple of smallest tstep')
# plt.ylabel('completeness')
# plot(tscales, mean_comp_walking, color=(0.,0.,1.), label='walking', lw=3.)
# plot(tscales, mean_comp_whole, color=(1.,0.,1.), label='whole', lw=3.)
# plot(tscales, mean_comp_restricted, color=(0.,1.,0.), label='restricted', lw=3.)
# plot(tscales, mean_comp_compensated, color=(1.,0.,0.), label='compensated', lw=3.)
# plot(tscales, completenesses_noE[0], color='k', label='no erosion', lw=3.)
# legend()
# figure(6)
# plt.gca().set_xscale('log')
# plt.xlabel('multiple of smallest tstep')
# plt.ylabel('completeness')
# for i in completenesses_compensated:
#     plot(tscales, i, color=gray)
# plot(tscales, mean_comp_compensated, 'k', lw=3.)
# title('compensated')
#
# figure(7)
# plt.gca().set_xscale('log')
# plt.xlabel('multiple of smallest tstep')
# plt.ylabel('completeness')
# for i in completenesses_walking:
#     plot(tscales, i, color=(0.35,0.35,1.), ls='-.', lw=1.)
# for i in completenesses_restricted:
#     plot(tscales, i, color=(0.45,1.,0.45), ls=':', lw=1.)
# for i in completenesses_compensated:
#     plot(tscales, i, color=(1.,0.55,0.55), ls='--', lw=1.)
# plot(tscales, mean_comp_walking, color=(0.,0.,1.), lw=3., label='walking')
# plot(tscales, mean_comp_whole, color=(1.,0.,1.), lw=3., label='whole')
# plot(tscales, mean_comp_restricted, color=(0.,1.,0.), lw=3., label='restricted')
# plot(tscales, mean_comp_compensated, color=(1.,0.,0.), lw=3., label='compensated')
# plot(tscales, completenesses_noE[0], color='k', lw=3., label='no erosion')
# legend()
