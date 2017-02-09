#!/bin/env python
import time
import os, sys, shutil, copy
import numpy as np
import pyflamestk.pareto as pareto

"""
Read a file and calculate the Pareto set.


"""
start_time = time.time()
fname_results_in = "20170201_results_10k.out"

sim_results = pareto.SimulationResults()
sim_results.read_simulation_results(fname_sims=fname_results_in)
sim_results.calculate_pareto_set()
sim_results.write_pareto_set()

print("{} second".format(time.time()-start_time))
