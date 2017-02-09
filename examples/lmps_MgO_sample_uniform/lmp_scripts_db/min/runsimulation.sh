#!/bin/bash
if [ -e "out.dat" ]
then
  rm out.dat
fi

# settings for hermes.mse.ufl.edu
# source /opt/intel/composer_xe_2011_sp1.8.273/bin/compilervars.sh 'intel64'
# export LAMMPS_BIN=/home/eugene/bin/lmp_intx_comb

# settings for blues.lcrc.anl.gov
export LAMMPS_BIN=/blues/gpfs/home/jlarson/software/lammps-17Nov16/src/lmp_serial

# run lammps simulation
$LAMMPS_BIN  -i in.min> out.dat 
