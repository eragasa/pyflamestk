#!/bin/env python
import time
import pyflamestk.pyposmat_mpi as pyposmat_mpi
"""

Authors:
Eugene J. Ragasa, eragasa@ufl.edu, 10/24/2016 - present

This script is a serial version of the ability of pyposmat to develop a 
Buckingham potential for magnesium oxide for MgO.

Required files:
    pyposmat.config
    pyposmat.potential
    pyposmat.qoi

Shell Variables:
    PYTHONPATH must include the pyflamestk module
    LAMMPS_BIN must provide the location of a serial version of lammps

"""
start_time = time.time()

# --- simulation parameters ---
n_simulations = 1000 # number of simulations
n_iterations = 10 # number of iterations

# supported cull_types are: percentile, pct_error
cull_type = 'percentile'  # cull by percentile
cull_param = 50.          # cull by parameter

# --- nothing in here should have to be changed

mcsampler = pyposmat_mpi.MpiIterativeSampler(n_iterations=n_iterations,
                                      n_simulations=n_simulations,
                                      cull_type=cull_type,
                                      cull_param = cull_param)
mcsampler.run()
print("{} second".format(time.time()-start_time))
