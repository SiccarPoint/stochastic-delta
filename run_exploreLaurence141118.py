from delta_obj2 import delta
import numpy as np
from matplotlib.pyplot import figure, plot, show, legend
import matplotlib.pyplot as plt
from copy import copy

# import the real params
# all are in standard units of meters; fluxes in m**3

iterations = 1  # 50

font = {'size': 16}

SL_trajectory = np.loadtxt('real_SL.txt')[1:]
Q_real = np.loadtxt('real_Qs.txt')[1:]
# use the 1st step as the initial condition...

mydelta = delta()
in_dict = mydelta.read_input_file('explore_inputs.txt')
nt = int(in_dict['nt'])

areas = (1., )
drifts = (0.01, )
completenesses = {1.: {}}
num_roll_pres = {1.: {}}

for area in areas:
    for drift in drifts:
        completeness_walk = []
        num_roll_pres_walk = []
        mydelta = delta()
        for param in ('activity_py', 'erosion_py_width', 'depo_py_width'):
            in_dict[param] = area
        in_dict['drift'] = drift

        for i in xrange(iterations):
            print 'iteration ', i
            tscales, completeness_walk = mydelta.execute(
                in_dict, SL_trajectory, completeness_records=completeness_walk,
                graphs=True, Q=Q_real, walking_erosion_depo=True)
            num_roll_pres_walk.append(mydelta.final_preserved.sum())
            # graphs=True will plot the strata, but don't turn it on unless
            # you are only doing a single run, not a param space exploration

        completenesses[area][drift] = copy(completeness_walk[0])
        num_roll_pres[area][drift] = copy(num_roll_pres_walk[0])
        show()

        times = mydelta.output_times
        traj = SL_trajectory
        preserved = mydelta.final_preserved
        rnode = mydelta.rnode
        strat_eta = mydelta.strat_eta
        shoreline = mydelta.shoreline_positions

        f = open("times1.txt", "w+")
        for i in range(nt):
            f.write('%s' % times[i])
            f.close
        g = open("traj1.txt", "w+")
        for i in range(nt):
            g.write('%s' % traj[i])
            g.close
        h = open("preserved1.txt", "w+")
        for i in range(nt):
            h.write('%s' % preserved[i])
            h.close
        j = open("rnode1.txt","w+")
        for i in range(nt):
             j.write('%s' % rnode[i])
             j.close
        k = open("strat_eta1.txt","w+")
        for i in range(nt):
             k.write('%s' % strat_eta[i])
             k.close
        l = open("shoreline1.txt","w+")
        for i in range(nt):
             l.write('%s' % shoreline[i])
             l.close





        # this dictionaries of dictionaries preserve the array of
        # completenesses at various timescales (get the correcponding array of
        # timescales with mydelta.full_completeness()[0]), and also an int
        # that represents the actual number of preserved rollovers for each run

        # You might need to think through whether this is done right for your
        # purposes.

        # Check the delta_obj2 code documentation to see the other data you
        # can pull from the model as it runs. Many will be useful (e.g., (x, z)
        # of rollovers).
