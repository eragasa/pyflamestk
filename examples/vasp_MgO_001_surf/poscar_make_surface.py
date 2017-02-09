import numpy as np
import pyflamestk.vasp as vasp

"""
Author: Eugene J. Ragasa
Makes a 3x3x3 supercell
"""


fname_poscar_in = 'MgO_unit.poscar'
unitcell = vasp.Poscar()
unitcell.read_file(fname_in=fname_poscar_in)

fname_poscar_out = 'MgO_sc112_s001.poscar'
sc_parameters = [1,1,3]
supercell = vasp.make_super_cell(unitcell,sc_parameters)

z_layers = []
n_atoms_old = 0
for z_max in np.linspace(0,1,101):
    n = 0
    for a in supercell.atoms:
        if a.position[2] < z_max: 
            n += 1
    n_atoms = n
    if n_atoms > n_atoms_old:
        z_layers.append(z_max)
    n_atoms_old = n_atoms
print("n_layers = {}".format(len(z_layers)))
print(z_layers)
supercell.write_file(fname_out=fname_poscar_out)


#fname_poscar_out = 'MgO_sc113_s001.poscar'
#sc_parameters = [1,1,3]
#supercell = vasp.make_super_cell(unitcell,sc_parameters)
#supercell.write_file(fname_out=fname_poscar_out)
