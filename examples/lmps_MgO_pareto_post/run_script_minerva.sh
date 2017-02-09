#!/bin/bash

# location of the pyflamestk module
export PYFLAMESTK_LIB=$(cd ../../;pwd)

# location of LAMMPS serial binary
export LAMMPS_BIN=/usr/local/bin/lammps

# build PYTHONPATH
# export PYTHONPATH=$PYFLAMESTK_LIB:$PYTHONPATH
export PYTHONPATH=$PYFLAMESTK_LIB

# requires python3 version of Anaconda from Continuum Analytics
python3 buckingham_pareto_report.py 
