# Usage: python abinit.py [trajectory] [run COLUMBUS?]
# trajectory number, to name the generated xyz geometry file for the deciding time point
# if you write anything after the trajectory number, the ab initio value will be calculated
#     otherwise, the already-calculated ab initio value will be used

# Import modules, read system argument variables, and define function
import math
import sys
import numpy as np
import os
import subprocess
print("\nModules imported\n")

if len(sys.argv) > 1:
    traj_num = int(sys.argv[1])
else:
    print("TRAJ NEEDED")
    exit()
if len(sys.argv) > 2:
    runcol = True
else:
    runcol = False

def compute_distance(coord1, coord2) -> float:
    '''Calculate distance between two Cartesian geometry points'''
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(coord1, coord2)))
stoptime = 0

# Check trajectory
data = []
direc = 'traj' + str(traj_num) + '/'
with open(direc + 'output.xyz', 'r') as file:
    lines = file.readlines()

i = 0
stop_processing = False

while i < len(lines):
    if lines[i].strip().isdigit():
        num_atoms = int(lines[i].strip())
        i += 1
        time_info = lines[i].strip().split()
        time = float(time_info[1])
        state1 = int(time_info[2])
        state2 = int(time_info[3])
        i += 1
        
        atom_data = {
            'time': time,
            'state1': state1,
            'state2': state2,
            'atoms': []
        }
        
        for _ in range(num_atoms):
            atom_line = lines[i].strip().split()
            atom = atom_line[0]
            coords = [float(coord) for coord in atom_line[1:]]
            atom_data['atoms'].append({'atom': atom, 'coords': coords})
            i += 1
        
        data.append(atom_data)
        
    else:
        i += 1

# Print the parsed data
for entry in data:

    # Read atom positions
    n_coords = None
    h_coords = []
    for atom_info in entry['atoms']:
        atom = atom_info['atom']
        coords = atom_info['coords']
        if atom == 'N':
            n_coords = coords
        elif atom == 'H':
            h_coords.append(coords)
    
    # Examine bond distances
    NHdists = [] # angstroms
    HHdists = []
    if n_coords is None or len(h_coords) < 3:
        print("Not enough atom coordinates.")
    else:
        # NH
        Hnum = 0
        for h_coord in h_coords:
            Hnum += 1
            distance = compute_distance(n_coords, h_coord)
            NHdists.append(distance)
            #print(f"Distance between N  and H{Hnum}: {distance:.6f} angstroms")
            # Terminate analysis
            if distance > 5.29:
                stop_processing = True
        # HH
        for h1 in range(3):
            for h2 in range(h1 + 1, 3):
                distance = compute_distance(h_coords[h1], h_coords[h2])
                HHdists.append(distance)
                #print(f"Distance between H{h1 + 1} and H{h2 + 1}: {distance:.6f} angstroms")
                # Terminate analysis
                if distance > 5.29:
                    stop_processing = True
        
        # Determine product channel
        channel = "UNDETERMINABLE" # needs to be manually examined
        if stop_processing: # If distance between any 2 atoms > 5.29 Å, TRAJECTORY IS OVER
            NHdists = sorted(NHdists)
            HHdists = sorted(HHdists)
            # If two N-H distances > 3 Å and if a N-H distance and H-H distance is under 2 Å
            if NHdists[1] > 3 and NHdists[0] < 2 and HHdists[0] < 2:
                channel = "molecular"
            # Elif two N-H distances < 2 Å
            elif NHdists[1] < 2:
                channel = "radical"
            # Save electronic state and quit
            print(f"Stopped due to distance > 10 a.u. at {entry['time']} fs")
            print(f"Channel: {channel}     diag: {entry['state1']}     MCH:  {entry['state2']}")
            print("NHdists (angstrom):   ", end="")
            for dist in NHdists:
                print(format(dist, ".2f") + "   ", end="")
            print("\nHHdists (angstrom):   ", end="")
            for dist in HHdists:
                print(format(dist, ".2f") + "   ", end="")
            print()
            # Print Hamiltonian (in a.u.)
            step = int(entry['time']*2)
            dat = open(direc + "output.dat", "r")
            datlines = dat.readlines()
            dat.close()
            for i in range(38 + step*62, 38 + step*62 + 5):
                for element in datlines[i][:-1].split():
                    print(format(float(element), "8.4f"), end="")
                print()
            stoptime = entry['time']
            print("Trajectory terminated at " + str(stoptime) + " fs")
            break

if not stop_processing: # no dissociation
    entry = data[-1]
    print(f"Stopped due to no dissociation by {entry['time']} fs")
    print(f"Result: NoDissociation     diag: {entry['state1']}     MCH:  {entry['state2']}")
    print("Trajectory did not terminate before " + str(stoptime) + " fs")

# Save geometry and model energy
def rip_from_traj(traj_num: int) -> dict:
    '''returns information from all time steps of a specified trajectory'''
    print("\nExamining trajectory " + str(traj_num))
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

traj_num = sys.argv[1]
# read output.lis
lis_data = rip_from_traj(traj_num)
# read output.dat
with open("traj" + str(traj_num)  + "/output.dat", "r") as file:
    lines = file.readlines()
time_limit = 1000 # fs
step = int(stoptime*2 + 1)
start_line = 35 + step*62
# Hamiltonian elements
H11 = float(lines[start_line + 3].split()[0]) # model S0 in Ha
H12 = float(lines[start_line + 3].split()[2]) 
H22 = float(lines[start_line + 4].split()[2]) # model S1 in Ha
H33 = float(lines[start_line + 5].split()[4])
# atom positions in Bohr radii
N  = [float(x) for x in lines[start_line + 53].split()]
H1 = [float(x) for x in lines[start_line + 54].split()]
H2 = [float(x) for x in lines[start_line + 55].split()]
H3 = [float(x) for x in lines[start_line + 56].split()]

# Run COLUMBUS
BUcolpath = '../../../COL/BU' # Path to the COLUMBUS calculation
colpath   = direc + 'COL'
if runcol:
    os.system("cp -r " + BUcolpath + " " + colpath)
    with open(colpath + '/geom', "w") as f:
        f.write(" N     7.0" + format(N[0],  "14.8f") + format(N[1],  "14.8f") + format(N[2],  "14.8f") + "   14.00307401\n")
        f.write(" H     1.0" + format(H1[0], "14.8f") + format(H1[1], "14.8f") + format(H1[2], "14.8f") + "    1.00782504\n")
        f.write(" H     1.0" + format(H2[0], "14.8f") + format(H2[1], "14.8f") + format(H2[2], "14.8f") + "    1.00782504\n")
        f.write(" H     1.0" + format(H3[0], "14.8f") + format(H3[1], "14.8f") + format(H3[2], "14.8f") + "    1.00782504\n")

    rv = subprocess.run(["sbatch", "script.sh"], cwd=colpath, capture_output=True)
    print("\n" + rv.stdout.decode('utf8'))
else:
    eadjust = 56.4759855158
    eV = 0.05512395
    au2cm = 219474.63
    with open(colpath + '/LISTINGS/ciudgsm.sp', "r") as f:
        lines = f.readlines()
    print(lines[-22])
    print(lines[-21])
    ab1 = float(lines[-22].split()[4]) # ab initio S0 in Ha, absolute
    ab2 = float(lines[-21].split()[4]) # ab initio S1 in Ha, absolute
    ab1 += eadjust # ab initio S0 in Ha, relative
    ab2 += eadjust # ab initio S1 in Ha, relative
    H11 -= 0.0551239574268122 # remove energy shift
    # in cm-1
    H11 *= au2cm
    H22 *= au2cm
    ab1 *= au2cm
    ab2 *= au2cm
    print("  ab initio      SURFGEN      diff")
    print(format(ab1, "10.0f"), format(H11, "10.0f"), format(H11-ab1, "10.0f"))
    print(format(ab2, "10.0f"), format(H22, "10.0f"), format(H22-ab2, "10.0f"))
