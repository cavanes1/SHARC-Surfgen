# Usage: python datarip.py [number of trajectories, default = 5000]
# Output: data.txt

# module import
import numpy as np
import math
from os import listdir
import os
import subprocess
import sys
import time
print("\nModules imported\n")

def rip_from_traj(traj_num: int) -> dict:
    '''returns information from all time steps of a specified trajectory'''
    sys.stdout.write("\rExamining trajectory " + str(curr))
    sys.stdout.flush()
    with open("traj" + str(curr)  + "/output.lis", "r") as file:
        lines = file.readlines()
    data = []
    step = 0
    for line in lines:
        if line[0:4] == "    ":
            if "******" in line:
                sys.exit("\nERROR on geometry (asterisks) of traj " + str(curr))
                break
            elif not int(line.split()[0]) == step:
                sys.exit("\nWARNING on time step " + str(step) + " of traj " + str(curr))
            data.append(line.split()[2:8])
            step += 1
    data = np.array(data)
    return {'diag': data[:, 0].astype(int), 
            'MCH':  data[:, 1].astype(int), 
            'kin':  data[:, 2].astype(float), 
            'pot':  data[:, 3].astype(float), 
            'tot':  data[:, 4].astype(float), 
            'amom': data[:, 5].astype(float)}

start_time = time.time()

# extract data from all trajectories
curr = 1
all_traj = {'diag': [], 'MCH': []} #, 'kin': [], 'pot': [], 'tot': [], 'amom': []}
max_traj = int(sys.argv[1]) if len(sys.argv) > 1 else 5000 #4987
while curr <= max_traj:
    this_traj = rip_from_traj(curr)
    for key in all_traj:
        all_traj[key].append(this_traj[key])
    curr += 1
for key in all_traj:
    all_traj[key] = np.array(all_traj[key])

# prepare data for population plot
nstates = 5
time_limit = 1000 # fs
fractions = {'diag': [[] for i in range(nstates)], 'MCH': [[] for i in range(nstates)]}
for step in range(time_limit*2 + 1):
    for representation in all_traj:
        time_slice = all_traj[representation][:, step]
        for i in range(nstates):
            fractions[representation][i].append(np.count_nonzero(time_slice == i + 1)/max_traj)

print("\n\nData extraction took %s seconds" % (time.time() - start_time))

# save to file
with open("data.txt", "w") as file:
    for representation in fractions:
        file.write(representation + '\n')
        for state in fractions[representation]:
            file.write(str(state)[1:-1] + '\n')
