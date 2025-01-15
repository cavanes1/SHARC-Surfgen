from collections import Counter

# read in output file
f = open("output.txt", "r")
lines = f.readlines()
f.close()

# identify branching ratio
results = {"radical": [[], []],
           "molecular": [[], []],
           "UNDETERMINABLE": [[], []],
           "NoDissociation": [[], []]}
for line in lines:
    if "MCH" in line:
        dummy, key, dummy, diag, dummy, MCH = line.split()
        results[key][0].append(int(diag))
        results[key][1].append(int(MCH))

# output result
for key in results:
    print("\n" + key + ": " + str(len(results[key][0])))
    diags = Counter(results[key][0])
    MCHs = Counter(results[key][1])
    print("        Diag   MCH")
    for state in range(1, 5 + 1):
        print("    " + str(state) + ":"  + format(diags[state], "6d") + format(MCHs[state], "6d"))
