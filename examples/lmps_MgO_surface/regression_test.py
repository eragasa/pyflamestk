#/bin/env python
import time
import pyflamestk.pyposmat as pyposmat
import pyflamestk.potentials as potentials

"""
This script creates a reference regression file for comparing the numerical results of pyflamestk
between two installations of lammps.

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

fname_ref = 'regress_minerva.out'            # references file

f = open(fname_ref)
lines = f.readlines()
f.close()

values_r = []
for i,line in enumerate(lines):
    if i == 0:
        names_r = [n.strip() for n in line.split(',')]
    elif i == 1:
        types_r = [t.strip() for t in line.split(',')]
    else:
        values_r.append([float(v) for v in line.split(',')])

p_index_1 = [i for i, x in enumerate(types_r) if x == 'param']
q_index_1 = [i for i, x in enumerate(types_r) if x == 'qoi']
e_index_1 = [i for i, x in enumerate(types_r) if x == 'err']


anls = pyposmat.AnlParameterSampler()
results_name   = None
results_type   = None
results_values = []

for value in values_r:
    param_dict = {}
    for idx in p_index_1:
        param_dict[names_r[idx]] = value[idx]

    value_name,value_type,value = anls.evaluate_one_parameter_set(param_dict)

    results_name = value_name
    results_type = value_type
    results_values.append(value)

p_index_2 = [i for i, x in enumerate(results_type) if x == 'param']
q_index_2 = [i for i, x in enumerate(results_type) if x == 'qoi']
e_index_2 = [i for i, x in enumerate(results_type) if x == 'err']

col_width = 15
print("-----PARAMETERS---------------")
print("".join([""] + [results_name[i].rjust(col_width) for i in p_index_2]))
for i,r in enumerate(results_values):
    #print("".join([i] + [str(r[i]).rjust(col_width) for i in p_index_2]))
    print("".join([str(i)] + [str(r[j]).rjust(col_width) for j in p_index_2]))

col_width = 20
print("-----MATERIAL PROPERTIES DIFFERENCES------")
print("".join([""] + [results_name[i].rjust(col_width) for i in q_index_2]))
for i in range(len(values_r)):
    diff = [results_values[i][j] - values_r[i][j] for j in q_index_2]
    print("".join([str(i)] + ["{:.2E}".format(d).rjust(col_width) for d in diff]))
col_width = 20
print("-----ERRORS DIFFERENCES-------------------")
print("".join([results_name[i].rjust(col_width) for i in e_index_2]))
for i,r in enumerate(results_values):
    diff = [results_values[i][j] - values_r[i][j] for j in e_index_2]
    print("".join([str(i)] + ["{:.2E}".format(d).rjust(col_width) for d in diff]))
