import copy
import pyflamestk.vasp as vasp
import pyflamestk.base as base
import pyflamestk.lammps as lammps

"""
Author: Eugene J. Ragasa
Makes a 3x3x3 supercell
"""
         
fname_poscar_in  = 'MgO_NaCl_001_s.vasp'
vasp_unit_cell = vasp.Poscar()
vasp_unit_cell.read_file(fname_in=fname_poscar_in)
vasp_unit_cell.normalize_h_matrix()

fname_lammps_out = 'MgO_NaCl_001_s.structure' 
lmps_unit_cell = lammps.StructureFile(obj=vasp_unit_cell)
lmps_unit_cell.write_file(fname_out=fname_lammps_out)

