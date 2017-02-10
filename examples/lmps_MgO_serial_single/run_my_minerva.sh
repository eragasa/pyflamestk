#!/bin/bash

PYTHONPATH=$(cd ../../;pwd):PYTHONPATH
LAMMPS_BIN=/usr/local/bin/lammps

python buckingham_single.py
