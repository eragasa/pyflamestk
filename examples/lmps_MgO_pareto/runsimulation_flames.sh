#!/bin/bash

export PYTHONPATH=$(cd ../../;pwd)
echo The PYTHONPATH is $PYTHONPATH

# this is required for FLAMES clusters
# LAMMPS_BIN must point to a valid serial version of LAMMPS
export LAMMPS_BIN=$(~/;pwd)/bin/lmp_serial

bash cleanup.sh
python buckingham_pareto.py 
