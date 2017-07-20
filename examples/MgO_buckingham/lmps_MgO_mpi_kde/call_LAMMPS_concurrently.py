import time
import pyflamestk.pyposmat_mpi as pyposmat_mpi

"""
This script runs a parameter sampling for the Buckingham potential using a
uniform distribution for each parameter in the parameter set and running
a lammps simulation.  This script does a uniform sampling of the parameter
ter space and uses MPI for parallelization.  Concurrency by providing a 
segregated file space for the workers.

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

# ---- begin script ----
random_seed = 0
n_simulations = 1000
fname_results_in = "20170201_culled_10k.out"
fname_results_out = 'results.out'

mpi_sampler = pyposmat_mpi.MpiPyPosmatEngine(is_restart = False,
                                             random_seed = random_seed)
mpi_sampler.sample_parameter_space(n_simulations = n_simulations,
                                   fname_results_out = fname_results_out,
                                   fname_results_in = fname_results_in,
                                   sampler_type = 'kde')

