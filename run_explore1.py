from delta_obj2 import delta
import numpy as np
from matplotlib.pyplot import figure, plot, show, legend
import matplotlib.pyplot as plt
from copy import copy

# import the real params
# all are in standard units of meters; fluxes in m**3

which_one = 'lo'
iterations = 1  # 50

font = {'size':16}
sect_comps = np.loadtxt('section_completenesses.txt')

SL_trajectory = np.loadtxt('real_SL.txt')[1:]
Q_real = np.loadtxt('real_Qs.txt')[1:]
# use the 1st step as the initial condition...

mydelta = delta()
in_dict = mydelta.read_input_file('explore_inputs.txt')
nt = int(in_dict['nt'])

areas = (0.4, )  # (0.025, 0.05, 0.1, 0.2, 0.4, 0.8)
lo_areas = (0.1, )  # (0.025, 0.05, 0.1)
hi_areas = (0.4, )  # (0.2, 0.4, 0.8)
drifts = (0.1, )  # (0.01, 0.02, 0.05, 0.1, 0.25, 0.5)
if which_one is 'lo':
    completenesses_loarea = {0.025: {}, 0.05: {}, 0.1: {}}
    num_roll_pres_loarea = {0.025: {}, 0.05: {}, 0.1: {}}
    which_areas = lo_areas
else:  # 'hi'
    completenesses_hiarea = {0.2: {}, 0.4: {}, 0.8: {}}
    num_roll_pres_hiarea = {0.2: {}, 0.4: {}, 0.8: {}}
    which_areas = hi_areas

for area in which_areas:
    for drift in drifts:
        completeness_walk = []
        num_roll_pres_walk = []
        mydelta = delta()
        for param in ('activity_py', 'erosion_py_width', 'depo_py_width'):
            in_dict[param] = area
        in_dict['drift'] = drift

        for i in xrange(iterations):
            print 'iteration ', i
            tscales, completeness_walk = mydelta.execute(in_dict, SL_trajectory,
            completeness_records=completeness_walk, graphs=True, Q=Q_real,      ##### set graphs=False when doing lots of runs
            walking_erosion_depo=True)
            num_roll_pres_walk.append(mydelta.final_preserved.sum())

        if which_one is 'lo':
            completenesses_loarea[area][drift] = copy(completeness_walk[0])
            num_roll_pres_loarea[area][drift] = copy(num_roll_pres_walk[0])
        else:
            completenesses_hiarea[area][drift] = copy(completeness_walk[0])
            num_roll_pres_hiarea[area][drift] = copy(num_roll_pres_walk[0])
