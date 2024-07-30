# Output: data.txt

# module import
import numpy as np
import math
from os import listdir
import os
import subprocess
import sys
print("\nModules imported\n")

# analyze present directories
path = "./"
dirs = [directory for directory in os.listdir(path) if os.path.isdir(path+directory)]

# extract energies
# total trajectory energy
curr = 1
go = True
diags = []
kins = []
pots = []
eners = []
while go:
    sys.stdout.write("\rExamining trajectory " + str(curr))
    sys.stdout.flush()
    f = open("traj" + str(curr)  + "/output.lis", "r")
    lines = f.readlines()
    f.close()
    dgst = []
    ener = []
    kin = []
    pot = []
    step = 0
    for line in lines:
        if line[0:4] == "    ":
            if not int(line.split()[0]) == step:
                print("WARNING on time step " + str(step) + " of traj " + str(curr))
            dgst.append(int(line.split()[2]))
            ener.append(float(line.split()[6]))
            kin.append(float(line.split()[4]))
            pot.append(float(line.split()[5]))
            step += 1
    diags.append(dgst)
    kins.append(kin)
    pots.append(pot)
    eners.append(ener)

    # exit loop
    curr += 1
    if "traj" + str(curr) not in dirs:
        go = False
        print("\nLast trajectory: " + str(curr - 1))

f = open("data.txt", "w")
for traj in diags:
    f.write(str(traj) + "\n")
for traj in pots:
    f.write(str(traj) + "\n")
f.close()
print("\n")
