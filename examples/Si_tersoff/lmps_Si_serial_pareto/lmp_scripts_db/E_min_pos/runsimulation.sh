#!/bin/bash
if [ -e "out.dat" ]
then
  rm out.dat
fi

THIS_HOSTNAME=$(hostname -f)
if [ "$THIS_HOSTNAME" = hermes.mse.ufl.edu ]; then
  # settings for hermes.mse.ufl.edu
  source /opt/intel/composer_xe_2011_sp1.8.273/bin/compilervars.sh 'intel64'
  export LAMMPS_BIN=/home/eugene/bin/lmp_intx_comb
elif [ "$THIS_HOSTNAME" = minerva ]; then
   # settings for Eugene's laptop
   export LAMMPS_BIN=/usr/local/bin/lammps
elif [ "$THIS_HOSTNAME" = blues2 ]; then
  # settings for blues.lcrc.anl.gov
  export LAMMPS_BIN=/blues/gpfs/home/jlarson/software/lammps-17Nov16/src/lmp_serial
fi

# run lammps simulation
${LAMMPS_BIN} -i in.min_pos > out.dat
