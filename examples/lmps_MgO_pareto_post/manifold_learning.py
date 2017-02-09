import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import NullFormatter
from sklearn import manifold, datasets

import pyflamestk.pareto as pareto

start_time = time.time()
fname_results_in = "20170201_results_10k.out"
fname_pareto_in = "20170201_pareto_10k.out"
is_load_previous_results = True
free_param_names = ['chrg_Mg', 'MgO_A', 'MgO_rho', 'OO_A', 'OO_rho', 'OO_C']
qoi_ref_values = {'MgO_NaCl.a0': 4.246,
                  'MgO_NaCl.c11': 277.00031,
                  'MgO_NaCl.c12': 91.67016,
                  'MgO_NaCl.c44': 144.00722,
                  'MgO_NaCl.B': 144.007722,
                  'MgO_NaCl.G': 153.4468767,
                  'MgO_NaCl.fr_a': 10.9781666,
                  'MgO_NaCl.fr_c': 8.98642095,
                  'MgO_NaCl.sch':5.067179685,
                  'MgO_NaCl.001s': 0.055950069}

#%% simulation results 

if is_load_previous_results is True:
    # load results
    sim_results = pareto.SimulationResults()
    sim_results.qoi_ref = qoi_ref_values 
    sim_results.read_simulation_results(fname_sims=fname_results_in,
                                        fname_pareto=fname_pareto_in)
    

else:
    # initialize simulation results
    sim_results = pareto.SimulationResults()
    sim_results.qoi_ref = qoi_ref_values 
    sim_results.read_simulation_results(fname_sims=fname_results_in)
    sim_results.calculate_pareto_set()
    sim_results.write_pareto_set()
    
#%% calculate culled pareto set
print(80*'-')
print('culled pareto set by percent error')
print(80*'-')
sim_results.calculate_culled_set('pct_error',100)
print(sim_results._perf_req)
print('n_results:',sim_results._results.shape[0])
print('n_pareto:',sim_results._pareto.shape[0])
print('n_culled:',sim_results._culled.shape[0])

#%% isomap

err_idx = [i for i,v in enumerate(sim_results.types) if v == 'err']
pareto_curve = sim_results._pareto[:,err_idx]
culled_pareto_curve = sim_results._culled[:,err_idx]

pareto_curve_t = manifold.Isomap(5,2).fit_transform(pareto_curve)
culled_pareto_curve_t = manifold.Isomap(5,2).fit_transform(culled_pareto_curve)

plt.figure(1)
plt.scatter(pareto_curve_t[:,0],
            pareto_curve_t[:,1], c='b', alpha=0.3, s=1)
plt.show()

plt.figure(2)
plt.scatter(culled_pareto_curve_t[:,0], 
            culled_pareto_curve_t[:,1], c='r', alpha=0.3, s=1)
plt.show()