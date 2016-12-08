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

n_simulations = 1000
start_time = time.time()
mc_sampler = pyposmat.MonteCarloParameterSampler()
mc_sampler.create_lammps_simulations()
mc_sampler.run(n_simulations = n_simulations)

print("{} second".format(time.time()-start_time))
