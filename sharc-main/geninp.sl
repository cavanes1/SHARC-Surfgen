#!/bin/bash -l

#SBATCH --partition=shared
#SBATCH -c 2
#SBATCH --time=24:00:00
#SBATCH --job-name=init

ml gcc/9.3.0
ml
python -u geninp.py > ../geninp.log
