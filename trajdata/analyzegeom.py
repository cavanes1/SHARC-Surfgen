import math

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
    #print(f"Time: {entry['time']}")
    #print(f"State 1: {entry['state1']}")
    #print(f"State 2: {entry['state2']}")
    #for atom_info in entry['atoms']:
    #    atom = atom_info['atom']
    #    coords = atom_info['coords']
    #    print(f"{atom}: {coords}")
    #print()

    n_coords = None
    h_coords = []
    for atom_info in entry['atoms']:
        atom = atom_info['atom']
        coords = atom_info['coords']
        if atom == 'N':
            n_coords = coords
        elif atom == 'H':
            h_coords.append(coords)
    
    if n_coords is None or len(h_coords) < 3:
        print("Not enough atom coordinates.")
    else:
        truncate = False
        for h_coord in h_coords:
            distance = compute_distance(n_coords, h_coord)
            #print(f"Distance between N and H: {distance:.6f}")
            if distance > 4:
                stop_processing = True
                break
        
        if stop_processing:
            print(entry['time'])
            print("Stopped due to distance > 4.")
            print(f"State 1: {entry['state1']}")
            print(f"State 2: {entry['state2']}")
            break
        
    #print()
