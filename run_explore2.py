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

areas = (0.025, 0.05, 0.1, 0.2, 0.4, 0.8)
drifts = (0.01, 0.02, 0.05, 0.1, 0.25, 0.5)
completenesses = {0.025: {}, 0.05: {}, 0.1: {}, 0.2: {}, 0.4: {}, 0.8: {}}
num_roll_pres = {0.025: {}, 0.05: {}, 0.1: {}, 0.2: {}, 0.4: {}, 0.8: {}}

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
                graphs=False, Q=Q_real, walking_erosion_depo=True)
            num_roll_pres_walk.append(mydelta.final_preserved.sum())
            # graphs=True will plot the strata, but don't turn it on unless
            # you are only doing a single run, not a param space exploration

        completenesses[area][drift] = copy(completeness_walk[0])
        num_roll_pres[area][drift] = copy(num_roll_pres_walk[0])
        # this dictionaries of dictionaries preserve the array of
        # completenesses at various timescales (get the correcponding array of
        # timescales with mydelta.full_completeness()[0]), and also an int
        # that represents the actual number of preserved rollovers for each run

        # You might need to think through whether this is done right for your
        # purposes.

        # Check the delta_obj2 code documentation to see the other data you
        # can pull from the model as it runs. Many will be useful (e.g., (x, z)
        # of rollovers).
