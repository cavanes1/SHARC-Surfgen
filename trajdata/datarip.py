# Usage: python datarip.py [number of trajectories, default = 5000]
# Output: data.txt, formatted with each line like 'diag', state 1, ..., state 5, 'MCH', state1, ..., state 5
#                   repeating for each channel

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

# read in output file
with open("output.txt", "r") as f:
    lines = f.readlines()
channel_results = []
trajcount = 0
for i, line in enumerate(lines):
    if " MCH:" in line:
        trajcount += 1
        dummy, channel, dummy, diag, dummy, MCH = line.split()
        diag, MCH, tm = int(diag), int(MCH), float(lines[i - 1].split()[-2])
        if channel == 'molecular':
            if MCH > 2:
                channel_results.append('molX')
            else:
                channel_results.append('mola')
        else:
            channel_results.append(channel)

# extract data from all trajectories
curr = 1
dict_list = {'all':            {'diag': [], 'MCH': []}, 
             'radical':        {'diag': [], 'MCH': []},
             'molX':           {'diag': [], 'MCH': []},
             'mola':           {'diag': [], 'MCH': []},
             'mol':            {'diag': [], 'MCH': []},
             'NoDissociation': {'diag': [], 'MCH': []},
             'UNDETERMINABLE': {'diag': [], 'MCH': []}}
max_traj = int(sys.argv[1]) if len(sys.argv) > 1 else trajcount
while curr <= max_traj:
    this_traj = rip_from_traj(curr)
    curr_channel = channel_results[curr - 1]
    for key in dict_list[curr_channel]:
        dict_list[curr_channel][key].append(this_traj[key]) # respective channel
        dict_list['all'][key].append(this_traj[key]) # all channels
        if curr_channel[:3] == 'mol':
            dict_list['mol'][key].append(this_traj[key]) # molecular channel (both)
    curr += 1
delete_list = []
for dictionary in dict_list:
    if dict_list[dictionary]['MCH'] == []:
        delete_list.append(dictionary)
    else:
        for key in dict_list[dictionary]:
            dict_list[dictionary][key] = np.array(dict_list[dictionary][key])
for dictionary in delete_list:
    del dict_list[dictionary]

# prepare data for population plot
nstates = 5
time_limit = 1000 # fs
for dictionary in dict_list:
    old = dict_list[dictionary]
    new = {'diag': [[] for i in range(nstates)], 'MCH': [[] for i in range(nstates)]}
    for step in range(time_limit*2 + 1):
        for representation in old:
            time_slice = old[representation][:, step]
            for i in range(nstates):
                new[representation][i].append(np.count_nonzero(time_slice == i + 1)/len(old['MCH']))
    dict_list[dictionary] = new

print("\n\nData extraction took %s seconds" % (time.time() - start_time))

# save to file
with open("data.txt", "w") as file:
    for dictionary in dict_list:
        for representation in dict_list[dictionary]:
            file.write(dictionary + ':   ' + representation + '\n')
            for state in dict_list[dictionary][representation]:
                file.write(str(state)[1:-1].replace(',', '') + '\n')
