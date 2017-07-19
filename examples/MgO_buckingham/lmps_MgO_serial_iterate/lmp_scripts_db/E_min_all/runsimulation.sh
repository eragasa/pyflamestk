#!/bin/bash
# This is an example bash runscript to run simulations within this
# directory.

THIS_HOSTNAME=$(hostname -f)
if [ "$THIS_HOSTNAME" = minerva ]; then
   # settings for Eugene's laptop
   export LAMMPS_BIN=/usr/local/bin/lammps
fi

# run lammps simulation
$LAMMPS_BIN  -i in.min> out.dat 
