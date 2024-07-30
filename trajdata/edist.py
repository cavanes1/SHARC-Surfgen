# Output: edist.png

# settings
mkplt = False # Make histogram
outinit = True # initconds file is outside the directories
deld = True # QM directories have been deleted

# module import
import numpy as np
import math
if mkplt:
    import matplotlib.pyplot as plt
from os import listdir
import os
import subprocess
import sys
print("\nModules imported")

# analyze present directories
path = "./"
dirs = [directory for directory in os.listdir(path) if os.path.isdir(path+directory)]

# extract energies
# total trajectory energy
curr = 1
go = True
kins = []
pots = []
eners = []
while go:
    sys.stdout.write("\rExamining trajectory " + str(curr))
    sys.stdout.flush()
    f = open("traj" + str(curr)  + "/output.lis", "r")
    lines = f.readlines()
    f.close()
    ener = float(lines[5].split()[6])
    kin = float(lines[5].split()[4])
    pot = float(lines[5].split()[5])
    kins.append(kin)
    pots.append(pot)
    eners.append(ener)

    # exit loop
    curr += 1
    if "traj" + str(curr) not in dirs:
        go = False

if not deld:
    # other energies
    Ekins = []
    Epots = []
    Etots = []
    inits = []
    conv = 27.211386246
    a0toA=0.52917721067
    print("Traj #      dEkin    dEpot    dEtot    dEsur")
    if outinit:
        f = open("../initconds", "r")
        lines = f.readlines()
        f.close()
    for i in range(curr - 1):
        currdir = "traj" + str(i + 1)
        if not outinit:
            f = open(currdir  + "/initconds", "r")
            lines = f.readlines()
            f.close()
        # read initial condition
        geom = []
        veloc = []
        for j in range(16+14*i+1,16+14*i+5):
            geom.append(lines[j][:60])
            veloc.append(lines[j][61:])
        # write QM.in
        f = open(currdir + "/QM/QM.in", "w")
        f.write('     4\n            1313004601')
        for j in range(len(geom)):
            f.write("\n" + geom[j][1])
            for val in geom[j][8:47].split():
                f.write(format(float(val)*a0toA, "13.7f"))
            f.write("   ")
            for val in veloc[j].split():
                f.write(format(float(val), "13.7f"))
        f.close()
        # run QM
        os.system("chmod +x " + currdir  + "/QM/runQM.sh")
        #os.system("cd " + currdir + ";QM/runQM.sh")
        rv = subprocess.run(["sh", "QM/runQM.sh"],cwd=currdir,capture_output=True)
        #h = open(str_target + '/runls', "w")
        #h.write(rv.stdout.decode('utf8'))
        #h.close()
        # read QM.out
        f = open(currdir  + "/QM/QM.out", "r")
        QMlines = f.readlines()
        f.close()
        QME = float(QMlines[3].split()[2])*conv
        # print energy data
        Ekin = float(lines[22+14*i].split()[1])*conv
        Epot = float(lines[24+14*i].split()[1])*conv
        Etot = float(lines[26+14*i].split()[1])*conv
        Ekins.append(Ekin)
        Epots.append(Epot)
        Etots.append(Etot)
        nstr = format(i + 1, "8.0f")
        kstr = format(Ekin - kins[i], "9.3f") #dEkin
        pstr = format(pots[i] - Epot - 5.8, "9.3f") #dEpot
        tstr = format(eners[i] - Etot - 5.8, "9.3f") #dEtot
        sstr = format(pots[i] - QME, "9.3f") #dEsur
        print(nstr + kstr + pstr + tstr + sstr)

# determine range
lowest = min(eners)
highest = max(eners)
print('\nRange: ' + str(lowest) + ' - ' + str(highest))
lo = math.floor(lowest)
hi = math.ceil(highest)
div = 0.25
binl = np.arange(lo, hi + div, div)


# make plot
if mkplt:
    #plt.figure(figsize=(10,5))
    plt.xlim(lo - 1, hi + 1)
    plt.xlabel('Total trajectory energy at t = 0 (eV)')
    plt.ylabel('Number of trajectories')
    plt.hist(eners, bins = binl)
    plt.savefig("edist.png")
    plt.close()
    print("\nPlot saved to edist.png")

print(eners)
