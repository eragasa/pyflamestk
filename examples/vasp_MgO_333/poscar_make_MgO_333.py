import copy
import pyflamestk.vasp as vasp
import pyflamestk.base as base

"""
Author: Eugene J. Ragasa
Makes a 3x3x3 supercell
"""
         
fname_poscar_in  = 'MgO_NaCl_unit.vasp'
unit_cell = vasp.Poscar()
unit_cell.read_file(fname_in=fname_poscar_in)
new_a1 = 5.0

# change the lattice parameter to 5
fname_poscar_out = 'MgO_NaCl_111.vasp' 
unit_cell.normalize_h_matrix()
old_a1 = unit_cell.lattice_parameter
unit_cell.lattice_parameter = new_a1
unit_cell.write_file(fname_out=fname_poscar_out)

fname_poscar_out = 'MgO_NaCl_333.vasp'
sc_p = [3,3,3]
supercell = vasp.make_super_cell(unit_cell,sc_p)
supercell.write_file(fname_out=fname_poscar_out)

# test_change_lattice_parameter
print('old_a1 = {}'.format(old_a1))
print('new_a1 = {}'.format(new_a1))

# supercell tests
assert sc_p[0], supercell.h_matrix[0,0]
assert sc_p[1], supercell.h_matrix[1,1]
assert sc_p[2], supercell.h_matrix[2,2]

# number of atoms test
assert unit_cell.n_atoms * sc_p[0] * sc_p[1] * sc_p[2], supercell.n_atoms
