import copy
import pyflamestk.vasp as vasp
import pyflamestk.base as base

"""
Author: Eugene J. Ragasa


"""

def make_defect_cell(unit_cell,sc_p,atoms_to_add,atoms_to_remove):
    pass    
fname_poscar_in  = 'MgO_NaCl_unit.vasp'
fname_poscar_out = 'MgO_NaCl_333_fr_a.vasp'
new_a1 = 5.0
sc_p = [3,3,3]
atoms_to_remove = [['O',[0.50,0.50,0.50]]]
atoms_to_add    = [['O',[0.75,0.75,0.75]]]

# change the lattice parameter value of the unit cell
unit_cell = vasp.Poscar()
unit_cell.read_file(fname_in=fname_poscar_in)
unit_cell.normalize_h_matrix()
unit_cell.lattice_parameter = new_a1

make_defect_cell(unit_cell,sc_p,atoms_to_add,atoms_to_remove)

super_cell = vasp.make_super_cell(unit_cell,sc_p)

print ('removing atoms')
for a in atoms_to_remove:
    symbol = a[0]
    position = a[1]
    print('\t',symbol,position)
    super_cell.remove_atom(symbol,position)
print ('adding atoms')
for a in atoms_to_add:
    symbol = a[0]
    position = a[1]
    print('\t',symbol,position)
    super_cell.add_atom(symbol,position)

n_atoms_added = len(atoms_to_add)
n_atoms_removed = len(atoms_to_remove)
n_atoms_perfect = unit_cell.n_atoms * sc_p[0] * sc_p[1] * sc_p[2]

print('uc:',unit_cell.lattice_parameter)
print('sc:',super_cell.lattice_parameter)
assert unit_cell.lattice_parameter, super_cell.lattice_parameter

for i in range(3):
    for j in range(3):
        if i == j:
            assert unit_cell.h_matrix[i,j] * sc_p[i], super_cell.h_matrix[i,j]
        else:
            print(i,j, super_cell.h_matrix[i,j])
assert n_atoms_perfect + n_atoms_added - n_atoms_removed, super_cell.n_atoms
super_cell.write_file(fname_out=fname_poscar_out)

