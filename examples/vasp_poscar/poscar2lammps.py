import pyflamestk.base as base
import pyflamestk.vasp as vasp
import pyflamestk.lammps as lammps

fname_poscar = 'POSCAR'
fname_lammps = 'structure.in'
sc_params = [3,3,3]
sym_order = ['Mg','O']
poscar = vasp.Poscar()
poscar.read_file(fname_poscar)

sc = base.make_super_cell(poscar,[3,3,3])

lmps = lammps.Structure(sc)
lmps.write(fname_lammps, sym_order,'charge')