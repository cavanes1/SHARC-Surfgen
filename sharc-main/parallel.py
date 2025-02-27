# To run: python parallel.py [first trajectory] [last trajectory]

# Module import
import os
import subprocess
import sys

# Bash script
text = '''#!/bin/bash -l

#SBATCH --partition=parallel
#SBATCH -c 2
#SBATCH --time=50:00:00
#SBATCH --job-name={ST}{excl}

ml intel/2020.1
echo "Will run trajectories {ST} to {EN}"
./run.sh {ST} {EN}'''

# Check running jobs
cmd = "squeue -u cavanes1"
rv = subprocess.run(cmd, cwd="./", shell=True, capture_output=True)
SLURMout = rv.stdout.decode()
print(SLURMout)
excludelist = set()
for line in SLURMout.split("\n"):
    if "c" in line:
        x = line.split()[-1]
        if x[0] == "c":
            excludelist.add(x)

excl_str = ""
count = 0
if not len(excludelist) == 0:
    excl_str += "\n#SBATCH -x "
    for node in excludelist:
        count += 1
        excl_str += node
        if not len(excludelist) == count:
            excl_str += ","

# Submit job
flname = "script" + sys.argv[1] + ".sh"
f = open(flname, "w")
f.write(text.format(excl=excl_str, ST=sys.argv[1], EN=sys.argv[2]))
f.close()
os.system("sbatch " + flname)
