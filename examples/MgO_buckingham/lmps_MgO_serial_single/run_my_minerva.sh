#!/bin/bash

export PYTHONPATH=$(cd ../../;pwd):PYTHONPATH
export LAMMPS_BIN=/usr/local/bin/lammps

python buckingham_single.py
