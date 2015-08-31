from delta_obj_2 import delta
import numpy as np

mydelta = delta()
ins = mydelta.read_input_file('test_inputs.txt')
nt = int(ins['nt'])
SL_trajectory = np.loadtxt('SL_z.txt')
completenesses = []
for i in xrange(10):
    tscales, completenesses = mydelta.execute('test_inputs.txt', SL_trajectory, completeness_records=completenesses, graphs=False, walking_erosion_depo=True)
figure(1)
plt.gca().set_xscale('log')
plt.xlabel('multiple of smallest tstep')
plt.ylabel('completeness')
for i in completenesses:
    plot(tscales, i)
show()
