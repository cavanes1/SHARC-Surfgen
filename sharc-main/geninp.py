# Before running, use wigner.py, then compile sharc-surfgen without field if targ is nonzero
# Run using geninp.sh (make sure to chmod +x sharc-surfgen.x)
# After running, compile sharc-surfgen with field if necessary, then submit SHARC with submit.sl or parallel.py

# settings
trajct = 5000 # number of trajectories to prepare
targ = 8.6 # desired total trajectory energy (eV), 0 if not relevant
tol = 0.1 # tolerance in total trajectory energy (eV)
veloc = "external" # can be "random X.X", or "external" if using Wigner distribution
randrot = True # randomize initial orientation

# module import
import random
import math
import os
import subprocess
import numpy as np
print("Modules imported")

# rotation, source: https://www.blopig.com/blog/2021/08/uniformly-sampled-3d-rotation-matrices/
def rotator(x):
    """Apply a random rotation in 3D, with a distribution uniform over the
    sphere.

    Arguments:
        x: vector or set of vectors with dimension (n, 3), where n is the
            number of vectors

    Returns:
        Array of shape (n, 3) containing the randomly rotated vectors of x,
        about the mean coordinate of x.

    Algorithm taken from "Fast Random Rotation Matrices" (James Avro, 1992):
    https://doi.org/10.1016/B978-0-08-050755-2.50034-8
    """

    def generate_random_z_axis_rotation():
        """Generate random rotation matrix about the z axis."""
        R = np.eye(3)
        x1 = np.random.rand()
        R[0, 0] = R[1, 1] = np.cos(2 * np.pi * x1)
        R[0, 1] = -np.sin(2 * np.pi * x1)
        R[1, 0] = np.sin(2 * np.pi * x1)
        return R

    # There are two random variables in [0, 1) here (naming is same as paper)
    x2 = 2 * np.pi * np.random.rand()
    x3 = np.random.rand()

    # Rotation of all points around x axis using matrix
    R = generate_random_z_axis_rotation()
    v = np.array([
        np.cos(x2) * np.sqrt(x3),
        np.sin(x2) * np.sqrt(x3),
        np.sqrt(1 - x3)
    ])
    H = np.eye(3) - (2 * np.outer(v, v))
    M = -(H @ R)
    x = x.reshape((-1, 3))
    mean_coord = np.mean(x, axis=0)
    return x @ M

# extract data from Wigner.py output
conv = 27.211386246
a0toA = 0.52917721067
geoms = []
velocs = []
if veloc == "external":
    f = open("../initconds", "r")
    lines = f.readlines()
    f.close()
    #os.system("mv initconds ..")
    Ninit = int(lines[1].split()[1])
    print("Ninit = " + str(Ninit))
    currinit = 0
    currtraj = 1
    while currtraj <= trajct: #Ninit
        currinit += 1
        if currinit > Ninit:
            exit("Not enough initial conditions")
        else:
            geom = []
            veloclst = []
            for j in range(14*currinit+3, 14*currinit+7):
                geom.append(lines[j][:60])
                veloclst.append(lines[j][61:])
            if randrot:
                vects = []
                for j in range(len(geom)):
                    vects.append([float(x) for x in geom[j][8:47].split()])
                    vects.append([float(x) for x in veloclst[j].split()])
                vects = rotator(np.array(vects))
            if not targ == 0: # if binning energies
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
                rv = subprocess.run(["sh", "QM/inQM.sh"],cwd="./",capture_output=True)
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
            for j in range(len(geom)):
                if randrot:
                    gls = ""
                    vls = ""
                    for k in vects[2*j]:
                        gls += format(k, "13.8f")
                    for k in vects[2*j+1]:
                        vls += format(k, "13.8f")
                    geoms[-1].append(geom[j][:8] + gls + geom[j][47:] + "\n")
                    velocs[-1].append(vls[1:] + "\n")
                else:
                    geoms[-1].append(geom[j] + "\n")
                    velocs[-1].append(veloclst[j])
        currtraj += 1

# create inputs
for i in range(1,trajct+1):
    u = random.random()
    rseed = math.floor(100000*u)
    inpstr = ''' geomfile    "geom"
 veloc       {veltxt}
 
 nstates      5
 actstates    5
 state        2 mch
 coeff        auto
 rngseed       {rngseed}
 
 ezero        0.0
 tmax         1000.0
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
