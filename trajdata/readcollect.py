# Outputs: readcollect.txt [number of trajectories, default is lesser of 5000 or number in output.txt]
#              formatted with each line like species, amount at each time step, all creation times (except for NH3), looping
#          terminal

# Module imports
from collections import Counter
import numpy as np
import sys

# read in output file
with open("output.txt", "r") as f:
    lines = f.readlines()

# identify branching ratio
results = {"radical":        [[], [], []],
           "molecular":      [[], [], []],
           "UNDETERMINABLE": [[], [], []],
           "NoDissociation": [[], [], []]}
NHX_times, NHa_times = [], []
max_traj = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
trajcount = 0
for i, line in enumerate(lines):
    if trajcount >= max_traj:
        break
    elif " MCH:" in line:
        trajcount += 1
        dummy, channel, dummy, diag, dummy, MCH = line.split()
        diag, MCH, time = int(diag), int(MCH), float(lines[i - 1].split()[-2])
        results[channel][0].append(diag)
        results[channel][1].append(MCH)
        results[channel][2].append(time)
        if channel == 'molecular':
            if MCH > 2:
                NHX_times.append(time)
            else:
                NHa_times.append(time)

# process data for species amount vs time
all_times = np.linspace(0, 1000, 2001)
products = []
NH3, NH2, NHX, NHa, FRG = trajcount, 0, 0, 0, 0
for time in all_times:
    FRG += results["UNDETERMINABLE"][2].count(time)
    NH2 += results["radical"][2].count(time)
    NHX += NHX_times.count(time)
    NHa += NHa_times.count(time)
    NH3 = trajcount - FRG - NH2 - NHX - NHa
    products.append([FRG, NH2, NHX, NHa, NH3])
products = np.array(products)/trajcount
product_evolution = {'FRG': [products[:, 0], results['UNDETERMINABLE'][2]],
                     'NH2': [products[:, 1], results['radical'][2]],
                     'NHX': [products[:, 2], NHX_times],
                     'NHa': [products[:, 3], NHa_times],
                     'NH3': [products[:, 4]]}
# save results to file
with open("readcollect.txt", "w") as file:
    for product in product_evolution:
        file.write(product + '\n')
        # species amount vs time
        file.write(str(list(product_evolution[product][0]))[1:-1].replace(',', '') + '\n')
        # times of creation of each species besides NH3
        # can be used to generate a plot of time for each product at each energy
        if len(product_evolution[product]) > 1:
            file.write(str(list(product_evolution[product][1]))[1:-1].replace(',', '') + '\n')

# output result
print("\nTotal: " + str(trajcount) + " trajectories")
denom = trajcount # might exclude some trajectories later
for channel in results:
    print("\n" + channel + ": " + str(len(results[channel][0])) + "  (" + str(len(results[channel][0])*100/trajcount) + "%)")
    diags = Counter(results[channel][0])
    MCHs  = Counter(results[channel][1])
    print("        Diag   MCH")
    for state in range(1, 5 + 1):
        print("    " + str(state) + ":"  + format(diags[state], "6d") + format(MCHs[state], "6d"))
print()
