# deletes trajectory directories that devolved into numerically extreme values

lastdir = 5010 # highest number of calculated trajectory directory

import subprocess
import os

rv = subprocess.check_output(['grep trajdata ../sharc-main/slurm-* | sort -u'], shell=True, text=True)
rv = rv.split('\n')
fails = []
for abc in rv[:-1]:
    fails.append(abc.split('/')[-4])

for d in fails:
    os.system(f"rm -r {d}")
    os.system(f"mv traj{lastdir} {d}")
    lastdir -= 1
