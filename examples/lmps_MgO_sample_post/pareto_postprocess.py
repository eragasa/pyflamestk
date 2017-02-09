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

sim_results = pyflamestk.pareto.SimulationResults()
sim_results.read_simulation_results(fname_in=fname_results_1)
sim_results.set_qoi_type(qoi_type)

#TODO: read this from configuration file
sim_results.qoi_values = {}
sim_results.qoi_values['MgO_NaCl_latt'] = 4.1212
sim_results.qoi_values['MgO_NaCl_B']    = 226.
sim_results.qoi_values['MgO_NaCl_G']    = 92.

sim_results.calculate_pareto_set(qoi_keys=qoi_keys)
sim_results.cull_by_percentile(pct_kept=pct_kept_each_iteration)
sim_results.calculate_parameter_estimates(param_list=param_list)
sim_results.calculate_qoi_estimates(qoi_keys=qoi_keys)

pyflamestk.pareto.resample_from_kernel_density(sim_results = sim_results, 
                                               param_list = param_list, 
                                               param_list_not_free = param_list_not_free)


#%% simulation results

#pyflamestk.pareto.make_histograms_parameters_combined(sim_results)
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
                        
                plt.text(sim_results_extern.np_all_sims[3,x_idx_extern],
                         sim_results_extern.np_all_sims[3,y_idx_extern],
                         'OURS',
                         bbox={'facecolor':'white', 'alpha':1.0, 'pad':2}) 
                         
                            
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

#fname_results_2 = "results_pub.out"          # original data                
#sim_results_extern = pyflamestk.pareto.SimulationResults()
#sim_results_extern.read_simulation_results(fname_in=fname_results_2)
#sim_results_extern.set_qoi_type(qoi_type)

#make_2d_pareto_plots_with_external_data(sim_results=sim_results,
#                                        sim_results_extern=sim_results_extern)

#------------------------------------------------------------------------------

