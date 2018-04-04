from delta_obj2 import delta
import numpy as np
from matplotlib.pyplot import figure, plot, show, legend
import matplotlib.pyplot as plt
from matplotlib import cm
from copy import deepcopy
import cPickle as pickle

# import the real params
# all are in standard units of meters; fluxes in m**3

num_pres_strata = 1600
channel_prop_median = 0.2  # arbitrary. Prob too high
font = {'size': 16, }
numnodes = 50000
delr = 0.05
ST = 0.01
SF = 0.65
theta = 1.571
drift = 0.1

dt_true = 200.
Uinit = 0.0005
flux_scaling = 0.5  # arbitrary accounting for basin size

# make a dummy counter for plotting
figcounter = 0

mydeltasims = []  # store each sim here

Q_options = [np.loadtxt('SPU5e-4Acc5dt200_12.30.44_sedfluxout.txt'),
             np.loadtxt('SPU5e-4Acc10dt200_12.30.44_sedfluxout.txt'),
             np.loadtxt('SPU5e-4Acc20dt200_12.30.44_sedfluxout.txt')]

# screwed up the normalization in the SP runs, so apply janky but robust
# normalization:
Q_options_sde = [np.loadtxt('U5e-4Acc5dt200_13.40.11_sedfluxout.txt'),
                 np.loadtxt('U5e-4Acc10dt200_13.40.11_sedfluxout.txt'),
                 np.loadtxt('U5e-4Acc20dt200_13.40.11_sedfluxout.txt')]

for SP, SDE in zip(Q_options, Q_options_sde):
    ratio = SDE[0]/SP[0]
    SP *= ratio

accel_options = [5., 10., 20.]


def load_mydelta_and_strat(fname, accel):
    if np.isclose(accel, 5.):
        Q_in = Q_options[0]
    elif np.isclose(accel, 10.):
        Q_in = Q_options[1]
    elif np.isclose(accel, 20.):
        Q_in = Q_options[2]
    else:
        raise NameError

    len_drift = 800.
    Q_real = np.concatenate((Q_in[0]*np.ones(400), Q_in,
                             Q_in[-1]*np.ones(len_drift)))
    Q_real *= flux_scaling
    # ^does this need adjusting for dt?

    # add the subsidence curve. Same idea
    SL_rate = np.concatenate((Uinit*np.ones(400), accel*Uinit*np.ones(1200)))
    SL_trajectory = np.cumsum(SL_rate)
    # DO NOT multiply by dt. Let's just scale everything to dt = 1 for now.
    # add an initial water depth if necessary here:
    SL_trajectory += 0.1

    # load it up
    f = open(fname, 'rb')
    mydelta = pickle.load(f)
    f.close()

    # do the plot
    color = [item/Q_real.max() for item in Q_real]
    for i in xrange(num_pres_strata):
        plot(mydelta.radial_distances, mydelta.strata[i, :],
             c=cm.plasma(color[i]))

    # pick the rollovers & add them to the fig:
    for i in xrange(num_pres_strata):
        # fit a (SF-ST) angle line to each point. rollover is point with
        # largest intersect
        # remove unmodified floor:
        notfloor = np.where(
            np.logical_not(np.isclose(mydelta.strata[i, :], 0.)))[0]
        c = ((SF - ST) * mydelta.radial_distances[notfloor] +
             mydelta.strata[i, notfloor])
        diff = np.diff(c)
        rollover_subindexes = np.where(
            np.diff((diff < 0.).astype(int)) == 1)[0]
        rollover_index = notfloor[rollover_subindexes]
        plot(mydelta.radial_distances[rollover_index],
             mydelta.strata[i, rollover_index], 'k,')


for (Q_in, accel) in zip(Q_options, accel_options):
    # add buffers before and after to get to nt = 1600
    len_drift = 800.
    drift_up = np.arange(len_drift, dtype=float)/len_drift*(
        accel*Uinit - Q_in[-1]) + Q_in[-1]

    Q_real = np.concatenate((Q_in[0]*np.ones(400), Q_in,
                             Q_in[-1]*np.ones(len_drift)))
    Q_real *= flux_scaling
    # ^does this need adjusting for dt?

    # add the subsidence curve. Same idea
    SL_rate = np.concatenate((Uinit*np.ones(400), accel*Uinit*np.ones(1200)))
    SL_trajectory = np.cumsum(SL_rate)
    # DO NOT multiply by dt. Let's just scale everything to dt = 1 for now.
    # add an initial water depth if necessary here:
    SL_trajectory += 0.1

    mydelta = delta()
    ins = {}
    # build the input dict:
    ins['nt'] = num_pres_strata
    ins['n'] = numnodes
    ins['delr'] = delr
    ins['delt'] = dt_true
    ins['Q'] = 0.  # superceded by Q input direct
    ins['ST'] = ST
    ins['SF'] = SF
    ins['theta'] = theta
    ins['activity_py'] = channel_prop_median
    ins['erosion_py_width'] = channel_prop_median
    ins['depo_py_width'] = channel_prop_median
    ins['drift'] = drift

    completenesses = []
    tscales, completenesses = mydelta.execute(
        ins, SL_trajectory, completeness_records=completenesses, graphs=False,
        Q=Q_real)

    figure(figcounter)
    color = [item/Q_real.max() for item in Q_real]
    for i in xrange(num_pres_strata):
        plot(mydelta.radial_distances, mydelta.strata[i, :],
             c=cm.plasma(color[i]))

    # pick the rollovers & add them to the fig:
    for i in xrange(num_pres_strata):
        # fit a (SF-ST) angle line to each point. rollover is point with
        # largest intersect
        # remove unmodified floor:
        notfloor = np.where(
            np.logical_not(np.isclose(mydelta.strata[i, :], 0.)))[0]
        c = ((SF - ST) * mydelta.radial_distances[notfloor] +
             mydelta.strata[i, notfloor])
        diff = np.diff(c)
        rollover_subindexes = np.where(
            np.diff((diff < 0.).astype(int)) == 1)[0]
        rollover_index = notfloor[rollover_subindexes]
        plot(mydelta.radial_distances[rollover_index],
             mydelta.strata[i, rollover_index], 'k,')

    f = open('mydeltaSP' + str(int(accel)) + '.save', 'wb')
    pickle.dump(mydelta, f, protocol=pickle.HIGHEST_PROTOCOL)
    f.close()
    figcounter += 1


# final_pres = mydelta.final_preserved
# completeness_subsampled = []
# for i in xrange(1):
#     condition = np.random.rand(final_pres.size) < (float(num_pres_strata)/nt)
#     new_pres_strata = np.logical_and(final_pres, condition)
#     tsc, comp = mydelta.full_completeness(record=new_pres_strata)
#     completeness_subsampled.append(comp.copy())

# figure(7)
# plt.gca().set_xscale('log')
# plt.xlabel('Timescale (s)')
# plt.ylabel('Completeness')
# for i in completeness_subsampled:
#     plot(tsc, i, '0.5', ls='--')
# plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
# plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
# plot([],[],'0.5', ls='--', label='resampled completenesses')
# plot([],[],'k', lw=3, label='real sections')
# plt.gca().set_ybound(0,1)
# plt.rc('font', **font)
# legend(loc=4)
#
# # now the restricted channel version
# completeness_surface_walk = []  # forced with the mean "deep channel" proportion
# completeness_synth_walk = []  # this one is forced with py (num_pres_strata-1)/nt == 0.0258
# completeness_surface_rand = []
# completeness_synth_rand = []
# num_roll_pres_surf_w = []
# num_roll_pres_synth_w = []
# num_roll_pres_surf_r = []
# num_roll_pres_synth_r = []
# synth_py = float((num_pres_strata-1)/nt)
# mydelta1 = delta()
# mydelta2 = delta()
# mydelta3 = delta()
# mydelta4 = delta()
# for i in xrange(30):
#     tscales, completeness_surface_walk = mydelta1.execute('real_inputs.txt', SL_trajectory,
#     completeness_records=completeness_surface_walk, graphs=False, Q=Q_real,
#     walking_erosion_depo=True)
#     num_roll_pres_surf_w.append(mydelta1.final_preserved.sum())
#     tscales, completeness_synth_walk = mydelta2.execute('real_inputs_frfixed.txt', SL_trajectory,
#     completeness_records=completeness_synth_walk, graphs=False, Q=Q_real,
#     walking_erosion_depo=True)
#     num_roll_pres_synth_w.append(mydelta2.final_preserved.sum())
#     tscales, completeness_surface_rand = mydelta3.execute('real_inputs.txt', SL_trajectory,
#     completeness_records=completeness_surface_rand, graphs=False, Q=Q_real,
#     restricted_channel_mass_conserved=True)
#     num_roll_pres_surf_r.append(mydelta3.final_preserved.sum())
#     tscales, completeness_synth_rand = mydelta4.execute('real_inputs_frfixed.txt', SL_trajectory,
#     completeness_records=completeness_synth_rand, graphs=False, Q=Q_real,
#     restricted_channel_mass_conserved=True)
#     num_roll_pres_synth_r.append(mydelta4.final_preserved.sum())
#     ### for some reason the initial topo breaks these!!! run it without...
#
# figure(8)
# for i in xrange(len(completeness_surface_walk)):
#     plot(tscales, completeness_surface_walk[i],'darkblue', ls='--')
#     #plot(tscales, completeness_synth_walk[i],'skyblue', ls='-')
#     plot(tscales, completeness_surface_rand[i],'firebrick', ls='--')
#     #plot(tscales, completeness_synth_rand[i],'lightsalmon', ls='-')
# plot([],[],'darkblue', ls='--', label='random walk')
# #plot([],[],'skyblue', ls='-', label='random walk, forced py')
# plot([],[],'firebrick', ls='--', label='no system memory')
# #plot([],[],'lightsalmon', ls='-', label='no system memory, forced py')
# plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
# plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
# plot([],[],'k', lw=3, label='real sections')
# plt.gca().set_xscale('log')
# plt.xlabel('Timescale (s)')
# plt.ylabel('Completeness')
# plt.rc('font', **font)
# legend(loc=4)
#
# figure('8b')
# for i in xrange(len(completeness_surface_walk)):
#     plot(tscales, completeness_surface_walk[i],'darkblue', ls='--')
#     plot(tscales, completeness_surface_rand[i],'firebrick', ls='--')
#     plot(tscales, completeness_synth_walk[i],'skyblue', ls='-')
#     #plot(tscales, completeness_synth_rand[i],'lightsalmon', ls='-')
# plot([],[],'darkblue', ls='--', label='random walk, real py')
# plot([],[],'firebrick', ls='--', label='no system memory, real py')
# plot([],[],'skyblue', ls='-', label='random walk, forced py')
# #plot([],[],'lightsalmon', ls='-', label='no system memory, forced py')
# plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
# plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
# plot([],[],'k', lw=3, label='real sections')
# plt.gca().set_xscale('log')
# plt.xlabel('Timescale (s)')
# plt.ylabel('Completeness')
# plt.rc('font', **font)
# legend(loc=4)
#
# figure(9)
# for i in xrange(len(completeness_surface_walk)):
#     plot(tscales, completeness_surface_walk[i],'darkblue', ls='--')
#     plot(tscales, completeness_synth_walk[i],'skyblue', ls='-')
#     plot(tscales, completeness_surface_rand[i],'firebrick', ls='--')
#     plot(tscales, completeness_synth_rand[i],'lightsalmon', ls='-')
# plot([],[],'darkblue', ls='--', label='random walk, real py')
# plot([],[],'skyblue', ls='-', label='random walk, forced py')
# plot([],[],'firebrick', ls='--', label='no system memory, real py')
# plot([],[],'lightsalmon', ls='-', label='no system memory, forced py')
# plot([],[],'k', lw=3, label='real sections')
# plt.gca().set_xscale('log')
# plt.xlabel('Timescale (s)')
# plt.ylabel('Completeness')
# legend(loc=4)
#
# figure(10)
# plt.gca().set_xscale('log')
# plt.xlabel('Timescale (s)')
# plt.ylabel('Completeness')
# for i in completenesses:
#     plot(tsc, i, '0.5', ls='--')
# plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
# plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
# plot([],[],'0.5', ls='--', label='simulated completeness')
# plot([],[],'k', lw=3, label='real sections')
# plt.gca().set_ybound(0,1)
# plt.rc('font', **font)
# legend(loc=4)
#
# figure(11)
# for i in xrange(mydelta.strat_eta.shape[0]):
#     plot(mydelta.radial_distances, mydelta.strat_eta[i,:], 'k')
# plt.xlabel('Radial distance')
# plt.ylabel('Height')
# plt.gca().set_xbound(0,3.5)
#
# figure(12)
# for i in xrange(len(completeness_surface_walk)):
#     plot(tscales, completeness_synth_walk[i],'skyblue', ls='-')
# plot([],[],'skyblue', ls='-', label='random walk, forced py')
# plot([],[],'k', lw=3, label='real sections')
# plot(sect_comps[54:,0], sect_comps[54:,1], 'k', lw=3)
# plot(sect_comps[54:,0], sect_comps[54:,2], 'k', lw=3)  # only >2000 s
# plt.gca().set_xscale('log')
# plt.xlabel('Timescale (s)')
# plt.ylabel('Completeness')
# legend(loc=4)
#
# figure(13)
# for i in xrange(mydelta2.strat_eta.shape[0]):
#     plot(mydelta2.radial_distances, mydelta2.strat_eta[i,:], 'k')
# plt.xlabel('Radial distance (m)')
# plt.ylabel('Height (m)')
# plt.gca().set_xbound(0,3.5)
#
# figure(14)
# mystrata = np.where(np.random.rand(mydelta.strat_eta.shape[0])
#                     <0.025831564048124558)[0]
# for i in mystrata:
#     plot(mydelta.radial_distances, mydelta.strat_eta[i,:], 'k')
# plt.xlabel('Radial distance')
# plt.ylabel('Height')
# plt.gca().set_xbound(0,3.5)
#
# figure(15)
# plot(sect_comps[54:,0], sect_comps[54:,1], 'k--', lw=3, label='Section A')
# plot(sect_comps[54:,0], sect_comps[54:,2], 'k:', lw=3, label='Section B')  # only >2000 s
# plt.gca().set_xscale('log')
# plt.xlabel('Timescale (s)')
# plt.ylabel('Completeness')
# plt.gca().set_ybound(0,1)
# plt.rc('font', **font)
# legend(loc=4)
#
# figure(16)
# plot(mydelta.output_times, SL_trajectory, 'k', lw=3)
# plt.xlabel('Time (s)')
# plt.ylabel('Water surface elevation (m)')
# plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
# plt.gca().set_ybound(0.,0.3)
# plt.rc('font', **font)
#
# show()
