#!/bin/env python
import pyflamestk.pyposmat as pyposmat
import pyflamestk.potentials as potentials
import pyflamestk.pareto as pareto
import time
import shutil
import os
import sys
"""
Eugene J. Ragasa, eragasa@ufl.edu, 10/24/2016

This script demonstrates the ability of pyposmat to develop a Buckingham potential for MgO.


"""
start_time = time.time()

# --- simulation parameters ---
n_simulations = 100
n_iterations = 10
pct_kept_each_iteration = 50

# --- file information ---
fname_results="results.out"
fname_params ="params.dat"

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

# some functions
def archive_information(idx_iterations,
                        fname_results_src="results.out",
                        fname_params_src="params.dat"):
    
    # create directory
    archive_dir_name      = "iter_{:03d}".format(idx_iterations)
    fname_results_dest    = "{}/sim_results_{:03d}.dat".format(archive_dir_name,
                                                               idx_iterations)
    fname_params_dest     = "{}/params_{:03d}.dat".format(archive_dir_name,
                                                          idx_iterations)

    if not os.path.exists(archive_dir_name):
        print("creating archive directory...")
        os.makedirs(archive_dir_name)
        shutil.copyfile(src = fname_results_src,
                        dst = fname_results_dest)
        shutil.copyfile(src = fname_params_src,
                        dst = fname_params_dest)

        print("\t directory name: {}".format(archive_dir_name))
        print("\t {} -> {}".format(fname_results_src,fname_results_dest))
        print("\t {} -> {}".format(fname_params_src,fname_params_dest))

# initialize
for idx_iterations in range(n_iterations):
    if not os.path.exists("iter_{:03d}".format(idx_iterations)):
        # generate simulations
        if idx_iterations == 0:
            mc_sampler = pyposmat.MonteCarloParameterSampler()
            mc_sampler.create_lammps_simulations()
            mc_sampler.run(n_simulations = n_simulations)
        else:
            file_sampler = pyposmat.FileParameterSampler()
            file_sampler.create_lammps_simulations()
            file_sampler.run(fname_params = fname_params,
                             fname_results = fname_results)

        # analyze results
        sim_results = pareto.SimulationResults()
        sim_results.read_simulation_results(fname_sims=fname_results)
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

        pareto.resample_from_kernel_density(sim_results = sim_results, 
                                            param_list = param_list, 
                                            param_list_not_free = param_list_not_free,
                                            fname_param = fname_params,
                                            n_resamples = n_simulations)

        # archive data
        archive_information(idx_iterations=idx_iterations)
    else:
        msg_template = "iteration {} has already been completed, skipping..."
        str_out = msg_template.format(idx_iterations)
        print(str_out)

print("{} second".format(time.time()-start_time))
