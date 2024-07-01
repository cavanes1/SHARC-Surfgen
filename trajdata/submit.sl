#!/bin/bash -l

#SBATCH --partition=shared
#SBATCH -c 1
#SBATCH --time=01:00:00
#SBATCH --job-name=sp

ml gcc/9.3.0
ml
python -u edist.py
