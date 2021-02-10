#!/bin/bash

#SBATCH -J lammpsExps               # Job name
#SBATCH -o job.%j.out         # Name of stdout output file (%j expands to jobId)
#SBATCH -N 5                  # Total number of nodes requested
#SBATCH -p regular             # Batch queue (use 'single' for one-node jobs)
#SBATCH -t 02:30:00           # Run time (hh:mm:ss) - 1.5 hours

# Launch MPI-based executable
# mpiexec -n {nproc} ./a.out
#lmp -in -partition 5x8 $HOME/gcmc/lammps/exp.mu.sic 
lmp -in $HOME/gcmc/lammps/exp.mu.sic
