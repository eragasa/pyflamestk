#!/bin/env python
import time
import pyflamestk.pyposmat as pyposmat
import pyflamestk.potentials as potentials

"""
This script runs a parameter sampling for the Buckingham potential using a
uniform distribution for each parameter in the parameter set and running
a lammps simulation

Required Directories:
  lmp_scripts_db  
  structure_db

Required Files:
  potential.mod
  pyposmat.config
  pyposmat.potential
  pyposmat.qoi
"""
param_dict = {}
param_dict['chrg_Mg'] = +2.0
param_dict['chrg_O']  = -2.0
param_dict['MgMg_A']   = 0.0 
param_dict['MgMg_rho'] = 0.5
param_dict['MgMg_C']   = 0.0
param_dict['MgO_A']    = 821.6
param_dict['MgO_rho']  = 0.3242
param_dict['MgO_C']    = 0.0
param_dict['OO_A']     = 2274.00 
param_dict['OO_rho']   = 0.1490
param_dict['OO_C']     = 27.88

n_simulations = 1000
start_time = time.time()
mc_sampler = pyposmat.PyPosmatEngine()
mc_sampler.sampler_type = 'uniform'
mc_sampler.sample_parameter_space(n_simulations)
#mc_sampler.evaluate_parameter_set(param_dict)
#mc_sampler.create_lammps_simulations()
#mc_sampler.run(n_simulations = n_simulations)

print("{} second".format(time.time()-start_time))
