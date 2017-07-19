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
cull_type = 'percentile'
cull_param = 50.

class IterativeSampler(pyposmat.PyPosmatEngine):

    def __init__(self,n_iterations,n_simulations,
                 cull_type,cull_param):
        pyposmat.PyPosmatEngine.__init__(self)
        self._n_iterations = n_iterations
        self._n_simulations = n_simulations
        self._iter_dir_format = "iter_{:03d}"
        self._fname_results_format = "results_{:03d}.out"
        self._fname_pareto_format = "pareto_{:03d}.out"
        self._fname_culled_format = "culled_{:03d}.out"
        self._cull_type = cull_type
        self._cull_param = cull_param

    @property
    def n_iterations(self):
        return self._n_iterations

    @n_iterations.setter
    def n_iterations(self, val):
        self._n_iterations = val

    @property
    def n_simulations(self):
        return self._n_simulations

    @n_simulations.setter
    def n_simulations(self, n_sims):
        self._n_simulations = n_sims

    def run(self,n_iterations=None):
        if n_iterations is not None:
            self._n_iterations = n_iterations

        for i_iter in range(self._n_iterations):
            self._log('starting iteration loop {}'.format(i_iter))
            fname_results_out = self._fname_results_format.format(i_iter)
            fname_pareto_out  = self._fname_pareto_format.format(i_iter)
            fname_culled_out  = self._fname_culled_format.format(i_iter)
            # generations
            if i_iter == 0:
                # use uniform sampling the first time around
                self.sample_parameter_space(n_simulations = self._n_simulations,
                                            fname_results = fname_results_out,
                                            sampler_type = 'uniform')
            else:
                # use the culled results from the previous iteration
                fname_results_in = self._fname_culled_format.format(i_iter)

                self.sample_parameter_space(n_simulations = self._n_simulations,
                                            fname_results = fname_results_out,
                                            sampler_type = 'kde',
                                            fname_results_in = fname_results_in)

            sim_results = pareto.SimulationResults()
            sim_results.read_simulation_results(fname_sims=fname_results_out)
            sim_results.calculate_pareto_set()
            sim_results.calculate_culled_set(self._cull_type,
                                             self._cull_param)
            sim_results.write_pareto_set(fname_pareto_out)
            sim_results.write_culled_set(fname_culled_out)

mcsampler = IterativeSampler(n_iterations=n_iterations,
                             n_simulations=n_simulations,
                             cull_type=cull_type,
                             cull_param = cull_param)
mcsampler.run()
print("{} second".format(time.time()-start_time))
