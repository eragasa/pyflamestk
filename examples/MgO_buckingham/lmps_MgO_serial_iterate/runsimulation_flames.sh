#!/bin/bash
export LAMMPS_BIN=$(cd ~/;pwd)/bin/lmp_serial
echo The LAMMPS_BIN is $LAMMPS_BIN
export PYTHONPATH=$(cd ../../../;pwd)
echo The PYTHONPATH is $PYTHONPATH

bash cleanup.sh
python buckingham_iterate.py
