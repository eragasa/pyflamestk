#!/bin/bash
export LAMMPS_BIN=$(~/;)/bin/lmp_serial
export PYTHONPATH=$(cd ../../../;pwd)
echo The PYTHONPATH is $PYTHONPATH

bash cleanup.sh
python buckingham_iterate.py
