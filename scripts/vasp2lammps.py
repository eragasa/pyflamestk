#!/usr/bin/env python
import pyflamestk.base
import pyflamestk.lammps
import pyflamestk.vasp
import argparse

parser = argparse.ArgumentParser(description="parses a POSCAR file and converts into a lammps structure file")
parser.add_argument('-i', action='store', 
                    dest='vasp_structure_file', 
                    help="input file, VASP POSCAR file")
parser.add_argument('-o', action='store', 
                    dest='lmps_structure_file',
                    help="output file, LAMMPS structure file")
parser.add_argument('-sc', action='store',
                    dest='sc_parameters',
                    help="make supercell, (e.g. 2x2x2)")
parser.add_argument('-sym', action='store',
                    dest='symbol_order',
                    help="order of atoms, (e.g. Zr,Si,O)")        


args = parser.parse_args()
vasp_structure_file = args.vasp_structure_file
lmps_structure_file = args.lmps_structure_file
sc_parameters = args.sc_parameters
symbol_order = args.symbol_order

print("vasp_structure_file={}".format(vasp_structure_file))
print("lammps_structure_file={}".format(lmps_structure_file))
print("supercell_parameters={}".format(sc_parameters))
print("symbol_order={}".format(symbol_order))

poscar = pyflamestk.vasp.Poscar()
poscar.load(vasp_structure_file)
unit_cell = poscar.structure

super_cell = []
if args.sc_parameters == '':
    pass
else:
    sc_parameters = [int(i) for i in args.sc_parameters.split('x')]
    super_cell = unit_cell.makeSupercell(sc_parameters)

lmps_structure = pyflamestk.lammps.StructureFile()
lmps_structure.structure = super_cell
lmps_structure.write(filename_out = lmps_structure_file,
                     symbol_list = symbol_order.split(","))
                     
