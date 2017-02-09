#!/bin/env python
import time
#import os, sys, shutil, copy
#import numpy as np
import pyflamestk.pareto as pareto
import pyflamestk.paretopost as paretopost
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas.tools.plotting import scatter_matrix

"""
Pareto post-processing report.  This has to be run on a laptop.  This should
be converted to a Jupyter notebook.
"""

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
    
#%% calculate Culled set

print(80*'-')
print('culled pareto set by percentile')
print(80*'-')
sim_results.calculate_culled_set('percentile',50)
sim_results.write_culled_set('culled.out')
print(sim_results._perf_req)

#%% calculate culled pareto set
print(80*'-')
print('culled pareto set by percent error')
print(80*'-')
sim_results.calculate_culled_set('pct_error',100)
print(sim_results._perf_req)
print('n_results:',sim_results._results.shape[0])
print('n_pareto:',sim_results._pareto.shape[0])
print('n_culled:',sim_results._culled.shape[0])

#%%

err_idx = [i for i,v in enumerate(sim_results.types) if v == 'err']
qoi_idx = [i for i,v in enumerate(sim_results.types) if v == 'qoi']
param_idx = [i for i,v, in enumerate(sim_results.types) if v == 'param']
free_param_idx = [sim_results.names.index(n) for n in free_param_names]
             
results_param_d = {sim_results.names[i]:sim_results.results[:,i] for i in param_idx}
results_param_df = pd.DataFrame(results_param_d)

results_qoi_d = {sim_results.names[i]:sim_results.results[:,i] for i in param_idx}
results_qoi_df = pd.DataFrame(results_qoi_d)

results_err_d = {sim_results.names[i]:sim_results.results[:,i] for i in err_idx}
results_err_df = pd.DataFrame(results_err_d)

results_abs_err_d = {sim_results.names[i]:np.abs(sim_results.results[:,i]) for i in err_idx}
results_abs_err_df = pd.DataFrame(results_abs_err_d)

pareto_param_d = {sim_results.names[i]:sim_results.pareto[:,i] for i in free_param_idx}
pareto_param_df = pd.DataFrame(pareto_param_d)

pareto_qoi_d = {sim_results.names[i]:sim_results.pareto[:,i] for i in param_idx}
pareto_qoi_df = pd.DataFrame(pareto_qoi_d)

pareto_err_d = {sim_results.names[i]:sim_results.pareto[:,i] for i in err_idx}
pareto_err_df = pd.DataFrame(pareto_err_d)

pareto_abs_err_d = {sim_results.names[i]:np.abs(sim_results.pareto[:,i]) for i in err_idx}
pareto_abs_err_df = pd.DataFrame(pareto_abs_err_d)

culled_param_d = {sim_results.names[i]:sim_results.culled[:,i] for i in free_param_idx}
culled_param_df = pd.DataFrame(culled_param_d)

culled_qoi_d = {sim_results.names[i]:sim_results.culled[:,i] for i in param_idx}
culled_qoi_df = pd.DataFrame(culled_qoi_d)

culled_err_d = {sim_results.names[i]:sim_results.culled[:,i] for i in err_idx}
culled_err_df = pd.DataFrame(culled_err_d)

culled_abs_err_d = {sim_results.names[i]:np.abs(sim_results.culled[:,i]) for i in err_idx}
culled_abs_err_df = pd.DataFrame(culled_abs_err_d)

scatter_matrix(pareto_abs_err_df,
               alpha=0.2, 
               figsize=(20,20),
               diagonal='kde')

scatter_matrix(culled_abs_err_df,
               alpha=0.2, 
               figsize=(20,20),
               diagonal='kde')

scatter_matrix(pareto_param_df,
               alpha=0.2, 
               figsize=(10,10),
               diagonal='kde')

scatter_matrix(culled_param_df,
               alpha=0.2, 
               figsize=(10,10),
               diagonal='kde')

#qoi_idx = [i for i,v in enumerate(sim_results.types) if v == 'qoi']

#results_err = results[:,[0] + err_idx]
#pareto_err = pareto[:,[0] + err_idx]
#results_abs_err = np.abs(results_err)
#pareto_abs_err = np.abs(pareto_err)

#results_err_df = pd.DataFrame({sim_results.names[i],sim_results.results[:,i] for i in err_idx})



#%% HISTOGRAMS OF ERRORS

#for i in err_idx:
#    plt.hist(sim_results.results[:,i])
#    plt.xlabel(sim_results.names[i])
#    plt.show()
    
#%% SCATTERPLOTS OF ERRORS
#for i in err_idx:
#    for j in err_idx:
#        if i < j:
#            plt.scatter(sim_results.results[i],
#                        sim_results.results[j])
#            plt.xlabel(sim_results.names[i])
#            plt.ylabel(sim_results.names[j])
#            plt.show()
            
#%% SCATTERPLOT OF ABS ERRORS
    
#%% HISTOGRAMS OF QOIS
#for i in qoi_idx:
    
#    plt.hist(sim_results.results[i])
#    plt.xlabel(sim_results.names[i])
#    plt.show()


#%% HISTOGRAMS OF PARAMETERS
#paretopost.simresults_create_param_histograms(sim_results,'results')
#paretopost.simresults_create_param_histograms(sim_results,'pareto')

#%% HISTOGRAMS OF QOIS
#print the values of a series
#qoi_name = 'MgO_NaCl.fr_a'
#data = sim_results.get_data(qoi_name,'pareto')
#print(np.isnan(data))
#print('qoi_name:',qoi_name,data.mean(),data.std())
#print(",\n".join([str(v) for v in sim_results.get_data('MgO_NaCl.fr_a','pareto')]))
#paretopost.simresults_create_qoi_histograms(sim_results,'pareto')

#%% HISTOGRAMS OF ERRORS
#paretopost.simresults_create_err_histograms(sim_results,'pareto')

#%%


#
print("{} second".format(time.time()-start_time))
