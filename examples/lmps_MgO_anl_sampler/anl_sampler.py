#!/bin/env python
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

# example parameter_set
param_names = ['chrg_Mg', 
           'OO_a',  'OO_c',  'OO_rho',
           'MgO_a', 'MgO_c', 'MgO_rho']
param_values = [2.0, -2.0, 0.0, 0.0, 0.5, 20764.00, 27.88, 0.1490,   821.6,  0.0, 0.3242]

param_dict = {}
for i in range(len(param_names)):
   param_dict[param_names[i]] = param_values[i] 

anls = pyposmat.PyPosmatEngine2()
value_name,value_type,value = anls.evaluate_parameter_set(param_dict)

print("p:",value_name)
print("q:",value_type)
print("e:",value)
