#!/bin/bash

export PYTHONPATH=$(cd ../../;pwd)
echo The PYTHONPATH is $PYTHONPATH

# tmux appears to be broken due to the brew install process on OSX
# this produces the error:
#     [warn] kq_init: detected broken kqueue; not using.: No such file or directory
# this comment suppresses the error
export EVENT_NOKQUEUE=1

export LAMMPS_BIN=/usr/local/bin/lammps
echo The LAMMPS_BIN is $LAMMPS_BIN

echo running cleanup script...
bash cleanup.sh
python buckingham_iterate.py
