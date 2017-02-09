import pyflamestk.vasp as vasp

fname_poscar_in = 'POSCAR'
fname_poscar_out = 'POSCAR.out'
sc_parameters = [3,3,3]
unitcell = vasp.Poscar()
unitcell.read_file(fname_in=fname_poscar_in)

supercell = vasp.make_super_cell(unitcell,sc_parameters)
supercell.write_file(fname_out=fname_poscar_out)

