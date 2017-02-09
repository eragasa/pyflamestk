import copy
import pyflamestk.vasp as vasp
import pyflamestk.base as base

"""
Author: Eugene J. Ragasa
Makes a 3x3x3 supercell
"""

fname_poscar_in  = 'MgO_NaCl_unit.vasp'
fname_poscar_out = 'MgO_NaCl_001_s.vasp'
min_slab_thickness = 20 # Angstroms
min_vac_thickness  = 20 #Angstroms
 
unit_cell = vasp.Poscar()
unit_cell.read_file(fname_in=fname_poscar_in)

unit_cell.normalize_h_matrix()

# minimum slab thickness
min_sim_box_z = min_slab_thickness + min_vac_thickness
# determine supercell size
a0 = unit_cell.lattice_parameter
h3 = unit_cell.h_matrix[2,:]
h3_norm = h3.dot(h3)**0.5
a3 = h3_norm * a0
min_sc_z = int(min_sim_box_z // a3 + 1) # floor division + 1
sc_p = [1,1,min_sc_z]
supercell = vasp.make_super_cell(unit_cell,sc_p)
# determine which atoms to remove
atoms_to_remove = []
for i,a in enumerate(supercell.atoms):
    if a.position[2] > (min_slab_thickness/supercell.a3) + 1e-6:
        atoms_to_remove.append(i)

# remove atoms in reverse order
for i in sorted(atoms_to_remove, reverse=True):
    del supercell.atoms[i]

supercell.write_file(fname_out=fname_poscar_out)
