# Usage: python onerip.py [trajectory number]
# Output: [trajectory number].txt, formatted with each line being a list of each of the following:
#             diag state, MCH state, E(S0), E(S1), E(T1), coupling, angmom, r(NH1), r(NH2), r(NH3)
#         terminal output describes when states change

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
    sys.stdout.write("\rExamining trajectory " + str(traj_num))
    sys.stdout.flush()
    with open("traj" + str(traj_num)  + "/output.lis", "r") as file:
        lines = file.readlines()
    data = []
    step = 0
    for line in lines:
        if line[0:4] == "    ":
            if "******" in line:
                sys.exit("\nERROR on geometry (asterisks) of traj " + str(traj_num))
                break
            elif not int(line.split()[0]) == step:
                sys.exit("\nWARNING on time step " + str(step) + " of traj " + str(traj_num))
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
traj_num = sys.argv[1]
# read output.lis
lis_data = rip_from_traj(traj_num)
# read output.dat
H11, H12, H22, H33, rNH1, rNH2, rNH3 = [], [], [], [], [], [], []
with open("traj" + str(traj_num)  + "/output.dat", "r") as file:
    lines = file.readlines()
time_limit = 1000 # fs
for step in range(time_limit*2 + 1):
    start_line = 35 + step*62
    H11.append(float(lines[start_line + 3].split()[0]))
    H12.append(float(lines[start_line + 3].split()[2]))
    H22.append(float(lines[start_line + 4].split()[2]))
    H33.append(float(lines[start_line + 5].split()[4]))
    Nx, Ny, Nz = lines[start_line + 53].split()
    H1x, H1y, H1z = lines[start_line + 54].split()
    H2x, H2y, H2z = lines[start_line + 55].split()
    H3x, H3y, H3z = lines[start_line + 56].split()
    rNH1.append(np.sqrt((float(Nx) - float(H1x))**2 + (float(Ny) - float(H1y))**2 + (float(Nz) - float(H1z))**2))
    rNH2.append(np.sqrt((float(Nx) - float(H2x))**2 + (float(Ny) - float(H2y))**2 + (float(Nz) - float(H2z))**2))
    rNH3.append(np.sqrt((float(Nx) - float(H3x))**2 + (float(Ny) - float(H3y))**2 + (float(Nz) - float(H3z))**2))
print("\n\nData extraction took %s seconds" % (time.time() - start_time))

# save to file
with open(f"{traj_num}.txt", "w") as file:
    file.write(str(list(lis_data['diag']))[1:-1].replace(',', '') + '\n')
    file.write(str(list(lis_data['MCH']))[1:-1].replace(',', '') + '\n')
    file.write(str(H11)[1:-1].replace(',', '') + '\n')
    file.write(str(H22)[1:-1].replace(',', '') + '\n')
    file.write(str(H33)[1:-1].replace(',', '') + '\n')
    file.write(str(H12)[1:-1].replace(',', '') + '\n')
    file.write(str(list(lis_data['amom']))[1:-1].replace(',', '') + '\n')
    file.write(str(rNH1)[1:-1].replace(',', '') + '\n')
    file.write(str(rNH2)[1:-1].replace(',', '') + '\n')
    file.write(str(rNH3)[1:-1].replace(',', '') + '\n')

# perform analysis and print to terminal
print("   ", lis_data['diag'][0], "   ", lis_data['MCH'][0], format(0, "8.1f"))
for step in range(1, time_limit*2 + 1):
    if (lis_data['diag'][step] != lis_data['diag'][step - 1]) or (lis_data['MCH'][step] != lis_data['MCH'][step - 1]):
        print("   ", lis_data['diag'][step], "   ", lis_data['MCH'][step], format(step/2, "8.1f"))
