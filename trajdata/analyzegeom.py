# Usage: python analyzegeom.py [trajectory] [verbosity]
verbosity = 2 # usually 1, except 2 when testing/debugging
traj = 0 # trajectory number, to name the generated xyz geometry file for the deciding time point
         # set to 0 if you do not want to generate this geometry file

import math
import sys

try:
    traj = int(sys.argv[1])
except:
    pass
try:
    verbosity = int(sys.argv[2])
except:
    pass

def compute_distance(coord1, coord2):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(coord1, coord2)))

data = []

with open('output.xyz', 'r') as file:
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
    # Print information for debugging
    if verbosity > 1:
        print(f"\nTime: {entry['time']} fs")
        print(f"diag: {entry['state1']}")
        print(f"MCH:  {entry['state2']}")
        for atom_info in entry['atoms']:
            atom = atom_info['atom']
            coords = atom_info['coords']
            print(f"{atom}: {coords}")
        print()

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
            if verbosity > 1:
                print(f"Distance between N  and H{Hnum}: {distance:.6f} angstroms")
            # Terminate analysis
            if distance > 5.29:
                stop_processing = True
        # HH
        for h1 in range(3):
            for h2 in range(h1 + 1, 3):
                distance = compute_distance(h_coords[h1], h_coords[h2])
                HHdists.append(distance)
                if verbosity > 1:
                    print(f"Distance between H{h1 + 1} and H{h2 + 1}: {distance:.6f} angstroms")
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
            if verbosity > 0:
                print("NHdists (angstrom):   ", end="")
                for dist in NHdists:
                    print(format(dist, ".2f") + "   ", end="")
                print("\nHHdists (angstrom):   ", end="")
                for dist in HHdists:
                    print(format(dist, ".2f") + "   ", end="")
                print()
                # Print Hamiltonian (in a.u.)
                step = int(entry['time']*2)
                dat = open("output.dat", "r")
                lines = dat.readlines()
                dat.close()
                #for i in range(35 + step*62, 35 + step*62 + 8):
                #    print(lines[i][:-1])
                for i in range(38 + step*62, 38 + step*62 + 5):
                    for element in lines[i][:-1].split():
                        print(format(float(element), "8.4f"), end="")
                    print()
            # Write geometry
            if traj > 0:
                step = int(entry['time']*2)
                f = open(str(traj) + ".xyz", "w")
                for i in range(step*6, step*6 + 6):
                    f.write(lines[i])
                f.close()
            break

if not stop_processing: # no dissociation
    entry = data[-1]
    print(f"Stopped due to no dissociation by {entry['time']} fs")
    print(f"Result: NoDissociation     diag: {entry['state1']}     MCH:  {entry['state2']}")
    # Write geometry
    if traj > 0:
        step = int(entry['time']*2)
        f = open(str(traj) + ".xyz", "w")
        for i in range(step*6, step*6 + 6):
            f.write(lines[i])
        f.close()
