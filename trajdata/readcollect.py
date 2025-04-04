from collections import Counter

# read in output file
with open("output.txt", "r") as f:
    lines = f.readlines()

# identify branching ratio
results = {"radical": [[], []],
           "molecular": [[], []],
           "UNDETERMINABLE": [[], []],
           "NoDissociation": [[], []]}
times, channels = [], []
trajcount = 0
for i, line in enumerate(lines):
    if " MCH:" in line:
        trajcount += 1
        dummy, channel, dummy, diag, dummy, MCH = line.split()
        results[channel][0].append(int(diag))
        results[channel][1].append(int(MCH))
        times.append(float(lines[i - 1].split()[-2]))
        channels.append(channel[0])

# output result
print("\nTotal: " + str(trajcount) + " trajectories")
denom = trajcount # might exclude some trajectories later
for channel in results:
    print("\n" + channel + ": " + str(len(results[channel][0])) + "  (" + str(len(results[channel][0])*100/trajcount) + "%)")
    diags = Counter(results[channel][0])
    MCHs = Counter(results[channel][1])
    print("        Diag   MCH")
    for state in range(1, 5 + 1):
        print("    " + str(state) + ":"  + format(diags[state], "6d") + format(MCHs[state], "6d"))
print()
