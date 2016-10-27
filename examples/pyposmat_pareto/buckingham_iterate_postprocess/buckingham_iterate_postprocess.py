#!/bin/env python
import pyflamestk.pyposmat
import pyflamestk.pareto
import matplotlib.pyplot as plt
import numpy as np

fname_results_1 = "results.out"      # data
fname_out       = "params.out"       # data

# --- parameter info ---
param_list = ['chrg_Mg',  'chrg_O',     
              'p_MgMg_a', 'p_MgMg_c', 'p_MgMg_rho', 
              'p_OO_a',   'p_OO_c',   'p_OO_rho',
              'p_MgO_a',  'p_MgO_c',  'p_MgO_rho']
# chrg_Mg = -chrg_O, no MgMg interaction
param_list_not_free = ['chrg_O', 
                       'p_MgMg_a', 'p_MgMg_c', 'p_MgMg_rho',
                       'p_MgO_c']

qoi_type      = "abserr"
qoi_keys     = ['MgO_NaCl_latt', 'MgO_NaCl_B', 'MgO_NaCl_G'] 
pct_kept_each_iteration = 50           

#TODO:  do this programmatically
n_iterations = 5

# do pareto analysis on iteration sets
sim_results = []
for idx in range(n_iterations):
    sim_results.append(pyflamestk.pareto.SimulationResults())
    fname_results = "iter_{:03d}/sim_results_{:03d}.dat".format(idx,idx)
    sim_results[idx].read_simulation_results(fname_in=fname_results)
    sim_results[idx].set_qoi_type(qoi_type)
    
    #TODO: read this from configuration file
    sim_results[idx].qoi_values = {}
    sim_results[idx].qoi_values['MgO_NaCl_latt'] = 4.1212
    sim_results[idx].qoi_values['MgO_NaCl_B']    = 226.
    sim_results[idx].qoi_values['MgO_NaCl_G']    = 92.
    
    sim_results[idx].calculate_pareto_set(qoi_keys=qoi_keys)
    sim_results[idx].cull_by_percentile(pct_kept=pct_kept_each_iteration)
    sim_results[idx].calculate_parameter_estimates(param_list=param_list)
    sim_results[idx].calculate_qoi_estimates(qoi_keys=qoi_keys)

import sys
sys.exit()

#%%
for i in range(n_iterations):
    for param in param_list:
        p_idx   = sim_results[i].all_names.index(param)
        p_value = np.mean(sim_results[i].np_all_sims[:,p_idx])
        p_std   = np.std(sim_results[i].np_all_sims[:,:p_idx])
        print("{} {} {} {}".format(i,param,p_value,p_std))

#%%
pyflamestk.pareto.make_histograms_parameters_combined(sim_results[0])
#%% simulation results_i

pyflamestk.pareto.make_histograms_parameters_combined(sim_results)
#pyflamestk.pareto.make_histograms_qois_combined(sim_results,qois=qoi_keys)
#pyflamestk.pareto.make_pareto_set_jpdf_plot(sim_results=sim_results,
#                                            param_list = param_list)
#pyflamestk.pareto.make_cull_set_jpdf_plot(sim_results=sim_results,
#                                          param_list = param_list)
#pyflamestk.pareto.make_simulations_jpdf_plot(sim_results=sim_results,
#                                             param_list = param_list)

# for making comparative 2D Pareto plots
# pyflamestk.pareto.make_2d_pareto_plots(sim_results=sim_results)


#%% pareto plots with external data
#TODO: need to generalize this
def make_2d_pareto_plots_with_external_data(sim_results,
                                            sim_results_extern,
                                            sim_results_iterate,
                                            i_iterate,
                                            labels_extern_short = "",
                                            labels_extern_long = "",
                                            show_dominated=True,
                                            show_pareto=True,
                                            show_cull=True):
    n = len(sim_results.qoi_keys)
    for i in range(n):
        for j in range(n):
            if i < j and i != j:
                x_label = sim_results.qoi_keys[i]
                y_label = sim_results.qoi_keys[j]
                print("{} {}".format(x_label,y_label))
    
                x_idx   = sim_results.all_names.index(x_label)
                y_idx   = sim_results.all_names.index(y_label)
                
                # apply performance requirements
                plt.figure()
                plt_handles = []
                if show_dominated:
                    plt_handles.append(plt.scatter(sim_results.np_all_sims[:,x_idx],
                                                   sim_results.np_all_sims[:,y_idx],
                                                   label='dominated',
                                                   marker='.',s = 1))
                if show_pareto:
                    plt_handles.append(plt.scatter(sim_results.np_pareto_set[:,x_idx],
                                                   sim_results.np_pareto_set[:,y_idx],
                                                   color='y',
                                                   label='Pareto'))
                if show_cull:
                    plt_handles.append(plt.scatter(sim_results.np_pareto_set_cull[:,x_idx],
                                                   sim_results.np_pareto_set_cull[:,y_idx],
                                                   color='r',
                                                   label='Pareto w/ constraints'))
                                                   
                if show_pareto:
                    pareto_front = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_all_sims[:,x_idx],
                                                                        sim_results.np_all_sims[:,y_idx],
                                                                        maxX = False, maxY= False)                
                    plt.plot(pareto_front[0],
                             pareto_front[1],
                             color='y',linewidth=2)
                             
                if show_cull:
                    pareto_front_cull = pyflamestk.pareto.pareto_frontier_2d(sim_results.np_pareto_set_cull[:,x_idx],
                                                                             sim_results.np_pareto_set_cull[:,y_idx],
                                                                             maxX = False, maxY= False)
                    plt.plot(pareto_front_cull[0],
                             pareto_front_cull[1],
                             color='r',linewidth=2)
                             
                # show external
                x_idx_extern   = sim_results_extern.all_names.index(x_label)
                y_idx_extern   = sim_results_extern.all_names.index(y_label)
                

                plt.text(sim_results_extern.np_all_sims[0,x_idx_extern],
                         sim_results_extern.np_all_sims[0,y_idx_extern],
                         'LC2.0',
                         bbox={'facecolor':'white', 'alpha':1.0, 'pad':2}) 
                        
                plt.text(sim_results_extern.np_all_sims[1,x_idx_extern],
                         sim_results_extern.np_all_sims[1,y_idx_extern],
                         'BG2.0',
                         bbox={'facecolor':'white', 'alpha':1.0, 'pad':2}) 
                        
                plt.text(sim_results_extern.np_all_sims[2,x_idx_extern],
                         sim_results_extern.np_all_sims[2,y_idx_extern],
                         'BG1.7',
                         bbox={'facecolor':'white', 'alpha':1.0, 'pad':2}) 
                        
                x_idx_iterate   = sim_results_iterate.all_names.index(x_label)
                y_idx_iterate   = sim_results_iterate.all_names.index(y_label)

                plt.text(sim_results_iterate.np_all_sims[i,x_idx_iterate],
                         sim_results_iterate.np_all_sims[i,y_idx_iterate],
                         "X",
                         bbox={'facecolor':'green', 'alpha':1.0, 'pad':2}) 
                         
                            
                plt.axis([0,
                          1.5*max(np.percentile(sim_results.np_pareto_set[:,x_idx],90),
                              max(sim_results_extern.np_all_sims[:,x_idx_extern])),
                          0,
                          1.5*max(np.percentile(sim_results.np_pareto_set[:,y_idx],90),
                               max(sim_results_extern.np_all_sims[:,y_idx_extern]))])
                plt.legend(handles=plt_handles)
                plt.xlabel(x_label)
                plt.ylabel(y_label)
                plt.show()         

# UNCOMMENT THIS OUT TO MAKE THIS WORK-------------------------------------------

fname_results_2 = "results_pub.out"          # original data                
sim_results_extern = pyflamestk.pareto.SimulationResults()
sim_results_extern.read_simulation_results(fname_in=fname_results_2)
sim_results_extern.set_qoi_type(qoi_type)

fname_results_3 = "results_iterate.out"          # original data                
sim_results_iterate = pyflamestk.pareto.SimulationResults()
sim_results_iterate.read_simulation_results(fname_in=fname_results_2)
sim_results_iterate.set_qoi_type(qoi_type)


make_2d_pareto_plots_with_external_data(sim_results=sim_results[2],
                                        sim_results_extern=sim_results_extern,
                                        sim_results_iterate=sim_results_iterate,
                                        i_iterate = 2)

#------------------------------------------------------------------------------

