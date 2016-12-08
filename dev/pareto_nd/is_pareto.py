# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 04:00:05 2016

@author: Eugene
"""

n_cycles = 10
fname_params      = "" 
fname_sim_results = ""
fname_log_file    = "log.pyposmat"
# first cycle

print("do initial simulation")
fname_sim_results = "results_000.dat"
print("\t  initial simulation file: {}".format(fname_sim_results))

for i_cycle in range(1,n_cycles):
    print("do pareto_analysis, iteration = {}".format(i_cycle))
    fname_params = "params_{:03d}.dat".format(i_cycle)
    print("\t parameter_file: {}".format(fname_params))