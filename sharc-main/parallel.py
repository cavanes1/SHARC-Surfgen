tot = 3000 # total number of trajectories
jobs = 10 # number of jobs to run

# Module import
import os

# Bash script
text = '''#!/bin/bash -l

#SBATCH --partition=parallel
#SBATCH -c 25
#SBATCH --time=36:00:00
#SBATCH --job-name=3RF{job}

ml intel/2020.1
./run.sh {ST} {EN} > run{job}.log'''

# Main code
for i in range(jobs):
    flname = "submit" + str(i + 1) + ".sl"
    f = open(flname, "w")
    size = tot//jobs
    f.write(text.format(job=str(i + 1), ST=str(i*size + 1), EN=str(i*size + size)))
    f.close()
    os.system("sbatch " + flname)
