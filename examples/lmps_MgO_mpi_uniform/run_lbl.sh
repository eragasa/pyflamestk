#!/bin/bash
rm slurm*
rm -rf mpi_worker_* 
sbatch -n 5 -p debug --time=30 lbl_concurrent_lammps.sl
