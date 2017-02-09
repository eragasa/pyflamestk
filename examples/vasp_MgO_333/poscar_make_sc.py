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
unit_cell.lattice_parameter = new_a1
unit_cell.write_file(fname_out=fname_poscar_out)

fname_poscar_out = 'MgO_NaCl_333.vasp'
sc_p = [3,3,3]
supercell = vasp.make_super_cell(unit_cell,sc_p)
supercell.write_file(fname_out=fname_poscar_out)


#------------------------------------------------------------------------------
# make 001 surface
#------------------------------------------------------------------------------
fname_poscar_out = 'MgO_NaCl_001_s.vasp'
min_slab_thickness = 10 # Angstroms
min_vac_thickness = 10 #Angstroms

# def make_001_surface(unit_cell, min_slab_thickness, min_vac_thickness):
assert issubclass(type(unit_cell), base.Structure)
# minimum slab thickness
min_sim_box_z = min_slab_thickness + min_vac_thickness
# determine supercell size
a0 = unit_cell.lattice_parameter
h3 = unit_cell.h_matrix[2,:]
h3_norm = h3.dot(h3)**0.5
a3 = h3_norm * a0
min_sc_z = min_sim_box_z // a3 + 1 # floor division + 1
sc_p = [1,1,min_sc_z]
supercell = vasp.make_super_cell(unit_cell,sc_p)

# determine which atoms to remove
atoms_to_remove = []
for i,a in enumerate(supercell.atoms):
    if a.position[2] > min_slab_thickness + 1e-6:
        atoms_to_remove.append(i)

# remove atoms in reverse order
for i in sorted(atoms_to_remove, reverse=True):
    del supercell.atoms[i]




#make_001_surface(unitcell, min_slab_thickness, min_vac_thickness)


supercell.write_file(fname_out=fname_poscar_out)

# return deep copy of supercell
#return copy.deepcopy(supercell)

#------------------------------------------------------------------------------
# make a Frenkel defect for an anion for MgO
#------------------------------------------------------------------------------

fname_poscar_out = 'MgO_NaCl_333_fr_a.vasp'
sc_p = [3,3,3]
atom_to_remove = ['O',[1,2,3]]
atom_to_add = ['O',[1,2,3]]
supercell = vasp.make_super_cell(unitcell,sc_p)
supercell.write_file(fname_out=fname_poscar_out)

#------------------------------------------------------------------------------
# make a Frenkel defect for a cation for MgO
#------------------------------------------------------------------------------

fname_poscar_out = 'MgO_NaCl_333_fr_c.vasp'
sc_p = [3,3,3]
atom_to_remove = ['Mg', [1,2,3]]
atom_to_add = ['Mg', [1,2,3]]
supercell = vasp.make_super_cell(unitcell,sc_p)
supercell.write_file(fname_out=fname_poscar_out)

#------------------------------------------------------------------------------
# make a Schottky defect for a cation for MgO
#------------------------------------------------------------------------------

fname_poscar_out = 'MgO_NaCl_333_sch.vasp'
sc_p = [3,3,3]
atom_to_remove = [['Mg',[1,2,3]],
                  ['O', [1,2,3]]]

supercell = vasp.make_super_cell(unitcell,sc_p)
supercell.write_file(fname_out=fname_poscar_out)

#------------------------------------------------------------------------------
# make 001 surface

