# Before running, use wigner.py, then compile sharc-surfgen without field
# Run using geninp.sh
# After running, compile sharc-surfgen with field if necessary, then submit SHARC with submit.sl or parallel.py

# settings
trajct = 3000 # number of trajectories to prepare
targ = 7.2 # desired total trajectory energy (eV), 0 if not relevant
tol = 0.1 # tolerance in total trajectory energy (eV)
veloc = "external" # can be "random X.X", or "external" if using Wigner distribution

# module import
import random
import math
import os
import subprocess
print("Modules imported")

# extract data from Wigner.py output
conv = 27.211386246
a0toA = 0.52917721067
geoms = []
velocs = []
if veloc == "external":
    f = open("initconds", "r")
    lines = f.readlines()
    f.close()
    os.system("mv initconds ..")
    Ninit = int(lines[1].split()[1])
    print("Ninit = " + str(Ninit))
    currinit = 0
    currtraj = 1
    while currtraj <= trajct: #Ninit
        currinit += 1
        if currinit > Ninit:
            exit("Not enough initial conditions")
        elif targ == 0:
            geoms.append([])
            velocs.append([])
            for j in range(14*currinit+3, 14*currinit+7):
                geoms[-1].append(lines[j][:60] + "\n")
                velocs[-1].append(lines[j][61:])
        else:
            geom = []
            veloclst = []
            for j in range(14*currinit+3,14*currinit+7):
                geom.append(lines[j][:60])
                veloclst.append(lines[j][61:])
            # write QM.in
            f = open("QM/QM.in", "w")
            f.write('     4\n            1313004601')
            for j in range(len(geom)):
                f.write("\n" + geom[j][1])
                for val in geom[j][8:47].split():
                    f.write(format(float(val)*a0toA, "13.7f"))
                f.write("   ")
                for val in veloclst[j].split():
                    f.write(format(float(val), "13.7f"))
            f.close()
            # run QM
            #os.system("chmod +x QM/runQM.sh")
            rv = subprocess.run(["sh", "QM/inQM.sh"],cwd="./",capture_output=True)
            #h = open(str_target + '/runls', "w")
            #h.write(rv.stdout.decode('utf8'))
            #h.close()
            # read QM.out
            f = open("QM/QM.out", "r")
            QMlines = f.readlines()
            f.close()
            QME = float(QMlines[3].split()[2])*conv
            # energy analysis
            Ekin = float(lines[8+14*currinit].split()[1])*conv
            TE = QME + Ekin
            print("Energy of initcond " + str(currinit) + " is " + format(TE, "9.3f"))
            if TE > targ + tol or TE < targ - tol:
                print("Energy out of range")
                continue
            print("Energy good, trajectory " + str(currtraj))
            geoms.append([])
            velocs.append([])
            for j in range(14*currinit+3, 14*currinit+7):
                geoms[-1].append(lines[j][:60] + "\n")
                velocs[-1].append(lines[j][61:])
        currtraj += 1

# create inputs
#f = open("inp/inp.all", "w")
for i in range(1,trajct+1):
    u = random.random()
    rseed = math.floor(100000*u)
    inpstr = ''' geomfile    "geom"
 veloc       {veltxt}
 
 nstates      2
 actstates    2
 state        2 mch
 coeff        auto
 rngseed       {rngseed}
 
 ezero        0.0
 tmax         500.0
 stepsize     0.5
 nsubsteps    25
 
 surf         sharc
 coupling     nacdr
 nogradcorrect
 ekincorrect  parallel_vel
 decoherence_scheme edc
 decoherence_param 0.2
 grad_all
 nac_select
 eselect      0.5\n'''.format(veltxt = veloc, rngseed = rseed)
    #f.write(inpstr)
    #g = open("inp/input." + str(i), "w")
    #g.write(inpstr)
    #g.close()
    os.system("mkdir -p ../trajdata/traj" + str(i))
    os.system("cp -r ./* ../trajdata/traj" + str(i))
    h = open("../trajdata/traj" + str(i) + "/input", "w")
    h.write(inpstr)
    h.close()
    if veloc == "external":
        g = open("../trajdata/traj" + str(i) + "/geom", "w")
        for line in geoms[i-1]:
            g.write(line)
        g.close()
        g = open("../trajdata/traj" + str(i) + "/veloc", "w")
        for line in velocs[i-1]:
            g.write(line)
        g.close()
    print("Writing trajectory directory " + str(i))
    rseed=rseed+1
#f.close()
