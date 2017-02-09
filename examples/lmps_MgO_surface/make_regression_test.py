#!/bin/env python
import time
import pyflamestk.pyposmat as pyposmat
import pyflamestk.potentials as potentials

"""
This script makes a regression test based on different parameter samplings
for a Buckingham potential.  Reference parmeter samplings.

Required Directories:
  lmp_scripts_db  
  structure_db

Required Files:
  potential.mod
  pyposmat.config
  pyposmat.potential
  pyposmat.qoi

References:
  [1] Lewis, G. V., and C. R. A. Catlow. Journal of Physics C: Solid State Physics 18.6 (1985): 1149.
  [2] Martinez, Jackelyn A., et al. Computer Physics Communications 203 (2016): 201-211.
  [3] Henkelman, Graeme, et al. Physical Review B 72.11 (2005): 115437.
"""


fname_regress_out = "regress_minerva.out"

# example parameter_set
param_names = ['chrg_Mg', 'chrg_O', 
           'p_MgMg_a', 'p_MgMg_c', 'p_MgMg_rho',
           'p_OO_a',  'p_OO_c',  'p_OO_rho',
           'p_MgO_a', 'p_MgO_c', 'p_MgO_rho']
pot_names = ['LC+2.0','BG+2.0','BG+1.7','MA+1.7']
potentials = {}
# Lewis Catlow, from Henkelman
potentials['LC+2.0'] = [2.0, -2.0, 
                          0.0, 0.0, 0.5, 
                          22764.00, 27.88, 0.1490,   
                          821.6,  0.0, 0.3242]
# Ball Grimes, from Henkelman
potentials['BG+2.0'] = [2.0, -2.0,
                          0.0, 0.0, 0.5,
                          9547.96, 32.0, 0.21916,
                          1279.69,  0.0, 0.29969]
# Ball Grimes, from Henkelman
potentials['BG+1.7'] = [1.7, -1.7,
                          0.0, 0.0, 0.5,
                          4870, 77.0, 0.2679,
                          929.69,0.0,0.29909]
# Martinez, parameterization A
potentials['MA+1.7'] = [1.7, -1.7, 
                          0.0, 0.0, 0.5,
                          22370,17.6997, 0.28323,
                          26007,299.981, 0.21938,]
# Martinez, parameterization B
#potentials['MB+1.7'] = [1.7, -1.7,
#                          0.0, 0.0, 0.5,
#                          33174,40.6299, 0.28429,
#                          70019,288.934, 0.29122]

anls = pyposmat.AnlParameterSampler()
results_name   = None
results_type   = None
results_values = []

for p in pot_names:
    param_dict = {}
    for i in range(len(param_names)):
        param_dict[param_names[i]] = potentials[p][i]

    value_name,value_type,value = anls.evaluate_one_parameter_set(param_dict)

    results_name = value_name
    results_type = value_type
    results_values.append(value)

p_index = [i for i, x in enumerate(results_type) if x == 'param']
q_index = [i for i, x in enumerate(results_type) if x == 'qoi']
e_index = [i for i, x in enumerate(results_type) if x == 'err']

col_width = 15
print("-----PARAMETERS---------------")
print("".join([""] + [results_name[i].rjust(col_width) for i in p_index]))
for i,r in enumerate(results_values):
    print("".join([pot_names[i]] + [str(r[i]).rjust(col_width) for i in p_index]))

col_width = 20
print("-----MATERIAL PROPERTIES------")
print("".join([""] + [results_name[i].rjust(col_width) for i in q_index]))
for i,r in enumerate(results_values):
    print("".join([pot_names[i]] + [str(r[i]).rjust(col_width) for i in q_index]))

col_width = 20
print("-----ERRORS-------------------")
print("".join([results_name[i].rjust(col_width) for i in e_index]))
for i,r in enumerate(results_values):
    print("".join([pot_names[i]] +[str(r[i]).rjust(col_width) for i in e_index]))

f = open(fname_regress_out,'w')
f.write(",".join(results_name)+"\n")
f.write(",".join(results_type)+"\n")
for r in results_values:
    # write results as csv line
    f.write(",".join([str(v) for v in r])+"\n")
f.close()
