#!/bin/bash -l

#SBATCH --partition=shared
#SBATCH -c 1
#SBATCH --time=01:00:00
#SBATCH --job-name=sp

ml intel/2020.1
./run.sh > run.log
