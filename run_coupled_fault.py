from delta_obj2 import delta
import numpy as np
from matplotlib.pyplot import figure, plot, show, legend
import matplotlib.pyplot as plt

# import the real params
# all are in standard units of meters; fluxes in m**3

num_pres_strata = 74
channel_prop_median = 0.18135607321131447  # using extract_ch_proportions.py
not_subaerial_prop_median = 0.61619800332778696
font = {'size':16}
sect_comps = np.loadtxt('section_completenesses.txt')

SL_trajectory = np.loadtxt('real_SL.txt')[1:]
Q_real = np.loadtxt('real_Qs.txt')[1:]
# use the 1st step as the initial condition...

mydelta = delta()
ins = mydelta.read_input_file('real_inputs.txt')
nt = int(ins['nt'])
completenesses = []
tscales, completenesses = mydelta.execute('real_inputs.txt', SL_trajectory,
completeness_records=completenesses, graphs=True, Q=Q_real,
initial_topo=(1.04, 0.065))  #, save_strat=(20,(3.5,0.4)))
final_pres = mydelta.final_preserved
completeness_subsampled = []
for i in xrange(10):
    condition = np.random.rand(final_pres.size)<(float(num_pres_strata-1)/nt)
    new_pres_strata = np.logical_and(final_pres, condition)
    tsc, comp = mydelta.full_completeness(record=new_pres_strata)
    completeness_subsampled.append(comp.copy())

figure(7)
plt.gca().set_xscale('log')
plt.xlabel('Timescale (s)')
plt.ylabel('Completeness')
for i in completeness_subsampled:
    plot(tsc, i, '0.5', ls='--')
plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
plot([],[],'0.5', ls='--', label='resampled completenesses')
plot([],[],'k', lw=3, label='real sections')
plt.gca().set_ybound(0,1)
plt.rc('font', **font)
legend(loc=4)

# now the restricted channel version
completeness_surface_walk = []  # forced with the mean "deep channel" proportion
completeness_synth_walk = []  # this one is forced with py (num_pres_strata-1)/nt == 0.0258
completeness_surface_rand = []
completeness_synth_rand = []
num_roll_pres_surf_w = []
num_roll_pres_synth_w = []
num_roll_pres_surf_r = []
num_roll_pres_synth_r = []
synth_py = float((num_pres_strata-1)/nt)
mydelta1 = delta()
mydelta2 = delta()
mydelta3 = delta()
mydelta4 = delta()
for i in xrange(30):
    tscales, completeness_surface_walk = mydelta1.execute('real_inputs.txt', SL_trajectory,
    completeness_records=completeness_surface_walk, graphs=False, Q=Q_real,
    walking_erosion_depo=True)
    num_roll_pres_surf_w.append(mydelta1.final_preserved.sum())
    tscales, completeness_synth_walk = mydelta2.execute('real_inputs_frfixed.txt', SL_trajectory,
    completeness_records=completeness_synth_walk, graphs=False, Q=Q_real,
    walking_erosion_depo=True)
    num_roll_pres_synth_w.append(mydelta2.final_preserved.sum())
    tscales, completeness_surface_rand = mydelta3.execute('real_inputs.txt', SL_trajectory,
    completeness_records=completeness_surface_rand, graphs=False, Q=Q_real,
    restricted_channel_mass_conserved=True)
    num_roll_pres_surf_r.append(mydelta3.final_preserved.sum())
    tscales, completeness_synth_rand = mydelta4.execute('real_inputs_frfixed.txt', SL_trajectory,
    completeness_records=completeness_synth_rand, graphs=False, Q=Q_real,
    restricted_channel_mass_conserved=True)
    num_roll_pres_synth_r.append(mydelta4.final_preserved.sum())
    ### for some reason the initial topo breaks these!!! run it without...

figure(8)
for i in xrange(len(completeness_surface_walk)):
    plot(tscales, completeness_surface_walk[i],'darkblue', ls='--')
    #plot(tscales, completeness_synth_walk[i],'skyblue', ls='-')
    plot(tscales, completeness_surface_rand[i],'firebrick', ls='--')
    #plot(tscales, completeness_synth_rand[i],'lightsalmon', ls='-')
plot([],[],'darkblue', ls='--', label='random walk')
#plot([],[],'skyblue', ls='-', label='random walk, forced py')
plot([],[],'firebrick', ls='--', label='no system memory')
#plot([],[],'lightsalmon', ls='-', label='no system memory, forced py')
plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
plot([],[],'k', lw=3, label='real sections')
plt.gca().set_xscale('log')
plt.xlabel('Timescale (s)')
plt.ylabel('Completeness')
plt.rc('font', **font)
legend(loc=4)

figure('8b')
for i in xrange(len(completeness_surface_walk)):
    plot(tscales, completeness_surface_walk[i],'darkblue', ls='--')
    plot(tscales, completeness_surface_rand[i],'firebrick', ls='--')
    plot(tscales, completeness_synth_walk[i],'skyblue', ls='-')
    #plot(tscales, completeness_synth_rand[i],'lightsalmon', ls='-')
plot([],[],'darkblue', ls='--', label='random walk, real py')
plot([],[],'firebrick', ls='--', label='no system memory, real py')
plot([],[],'skyblue', ls='-', label='random walk, forced py')
#plot([],[],'lightsalmon', ls='-', label='no system memory, forced py')
plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
plot([],[],'k', lw=3, label='real sections')
plt.gca().set_xscale('log')
plt.xlabel('Timescale (s)')
plt.ylabel('Completeness')
plt.rc('font', **font)
legend(loc=4)

figure(9)
for i in xrange(len(completeness_surface_walk)):
    plot(tscales, completeness_surface_walk[i],'darkblue', ls='--')
    plot(tscales, completeness_synth_walk[i],'skyblue', ls='-')
    plot(tscales, completeness_surface_rand[i],'firebrick', ls='--')
    plot(tscales, completeness_synth_rand[i],'lightsalmon', ls='-')
plot([],[],'darkblue', ls='--', label='random walk, real py')
plot([],[],'skyblue', ls='-', label='random walk, forced py')
plot([],[],'firebrick', ls='--', label='no system memory, real py')
plot([],[],'lightsalmon', ls='-', label='no system memory, forced py')
plot([],[],'k', lw=3, label='real sections')
plt.gca().set_xscale('log')
plt.xlabel('Timescale (s)')
plt.ylabel('Completeness')
legend(loc=4)

figure(10)
plt.gca().set_xscale('log')
plt.xlabel('Timescale (s)')
plt.ylabel('Completeness')
for i in completenesses:
    plot(tsc, i, '0.5', ls='--')
plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
plot([],[],'0.5', ls='--', label='simulated completeness')
plot([],[],'k', lw=3, label='real sections')
plt.gca().set_ybound(0,1)
plt.rc('font', **font)
legend(loc=4)

figure(11)
for i in xrange(mydelta.strat_eta.shape[0]):
    plot(mydelta.radial_distances, mydelta.strat_eta[i,:], 'k')
plt.xlabel('Radial distance')
plt.ylabel('Height')
plt.gca().set_xbound(0,3.5)

figure(12)
for i in xrange(len(completeness_surface_walk)):
    plot(tscales, completeness_synth_walk[i],'skyblue', ls='-')
plot([],[],'skyblue', ls='-', label='random walk, forced py')
plot([],[],'k', lw=3, label='real sections')
plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
plt.gca().set_xscale('log')
plt.xlabel('Timescale (s)')
plt.ylabel('Completeness')
legend(loc=4)

figure(13)
for i in xrange(mydelta2.strat_eta.shape[0]):
    plot(mydelta2.radial_distances, mydelta2.strat_eta[i,:], 'k')
plt.xlabel('Radial distance (m)')
plt.ylabel('Height (m)')
plt.gca().set_xbound(0,3.5)

figure(14)
mystrata = np.where(np.random.rand(mydelta.strat_eta.shape[0])
                    <0.025831564048124558)[0]
for i in mystrata:
    plot(mydelta.radial_distances, mydelta.strat_eta[i,:], 'k')
plt.xlabel('Radial distance')
plt.ylabel('Height')
plt.gca().set_xbound(0,3.5)

figure(15)
plot(sect_comps[54:,0], sect_comps[54:,1], 'k--', lw=3, label='Section A')
plot(sect_comps[54:,0], sect_comps[54:,2], 'k:', lw=3, label='Section B')  # only >2000 s
plt.gca().set_xscale('log')
plt.xlabel('Timescale (s)')
plt.ylabel('Completeness')
plt.gca().set_ybound(0,1)
plt.rc('font', **font)
legend(loc=4)

figure(16)
plot(mydelta.output_times, SL_trajectory, 'k', lw=3)
plt.xlabel('Time (s)')
plt.ylabel('Water surface elevation (m)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.gca().set_ybound(0.,0.3)
plt.rc('font', **font)

show()
