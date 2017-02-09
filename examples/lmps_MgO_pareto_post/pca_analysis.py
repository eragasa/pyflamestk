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

#%% pca analysis

n_rows, n_columns = sim_results.results.shape
param_idx = [i for i,v in enumerate(sim_results.types) if v == 'param']
qoi_idx = [i for i,v in enumerate(sim_results.types) if v == 'qoi']
err_idx = [i for i,v in enumerate(sim_results.types) if v == 'err']

param_names = [sim_results.names[i] for i in param_idx]
qoi_names = [sim_results.names[i] for i in qoi_idx]           
err_names = [sim_results.names[i] for i in err_idx]

culled_err = sim_results.culled[:,err_idx]
n_rows, n_cols = culled_err.shape
culled_err_means = np.array([np.mean(sim_results.culled[:,i]) for i in err_idx])
culled_err_scatter_matrix = np.zeros((n_cols,n_cols))
for i in range(n_rows):
    diff_v = culled_err[i,:].T - culled_err_means
    diff_v = diff_v.reshape((n_cols,1))
    culled_err_scatter_matrix += (diff_v).dot(diff_v.T)
print('culled_err_scatter_matrix:',culled_err_scatter_matrix.shape,'\n',
      culled_err_scatter_matrix)

culled_err_cov_matrix = np.cov(culled_err.T)

eig_val_cov, eig_vec_cov = np.linalg.eig(culled_err_cov_matrix)
eig_pairs = [(np.abs(eig_val_cov[i]), eig_vec_cov[:,i]) for i in range(len(eig_val_cov))]
eig_pairs.sort(key = lambda x: x[0], reverse = True)
for eig_pair in eig_pairs:
    print('eig_value:',eig_pair[0])
    print('eig_pair:',eig_pair[1])