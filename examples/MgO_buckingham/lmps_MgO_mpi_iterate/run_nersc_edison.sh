#!/bin/bash
rm slurm*
rm -rf mpi_worker_* 
sbatch -N 1 -p debug --time=30 runjob_nersc_edison.sl
