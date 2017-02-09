#!/bin/env python
import time
import pyflamestk.pyposmat as pyposmat
import pyflamestk.potentials as potentials

"""
This script runs a parameter sampling for the Buckingham potential using a
uniform distribution for each parameter in the parameter set and running
a lammps simulation.

Required Directories:
  lmp_scripts_db  
  structure_db

Required Files:
  potential.mod
  pyposmat.config
  pyposmat.potential
  pyposmat.qoi

Output Files:
  pyposmat.log - log file
  results.out - results file

Modify for the system:
  1.  lmp_scripts_db/*/runsimulation.sh provides the location of the lammps
  binary.  Alternatively, $LAMMPS_BIN variable can be set.
  2.  the pyflamestk module needs to be added to the $PYTHONPATH.
      export PYTHONPATH=$(cd ..;pwd):$PYTHONPATH
"""

n_simulations = 100
n_seed = None
is_restart = True
start_time = time.time()
mc_sampler = pyposmat.PyPosmatEngine2(is_restart = True, random_seed = None)
mc_sampler.sampler_type = 'uniform'
mc_sampler.sample_parameter_space(n_simulations)

print("{} second".format(time.time()-start_time))
