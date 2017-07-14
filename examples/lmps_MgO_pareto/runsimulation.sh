#!/bin/bash

export PYTHONPATH=$(cd ../../;pwd)
echo The PYTHONPATH is $PYTHONPATH
python buckingham_pareto.py 
