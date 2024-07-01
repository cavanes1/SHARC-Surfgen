# read in output file
f = open("output.txt", "r")
lines = f.readlines()
f.close()

# identify branching ratio
s11 = 0
s12 = 0
s21 = 0
s22 = 0
fails = 0
for i in range(len(lines)):
    if lines[i] == "State 1: 1\n":
        s11 += 1
    elif lines[i] == "State 1: 2\n":
        s12 += 1
    elif lines[i] == "State 2: 1\n":
        s21 += 1
    elif lines[i] == "State 2: 2\n":
        s22 += 1
    elif "entering" in lines[i]:
        if "entering" in lines[i+1]:
            fails += 1

# output result
print("\nState 1: 1 x " + str(s11))
print("State 1: 2 x " + str(s12))
print("State 2: 1 x " + str(s21))
print("State 2: 2 x " + str(s22))
print("Did not dissociate: " + str(fails))
print("\nDiagonal basis")
print("     State 1: " + str(s11) + "       State 2: " + str(s12) + "       No dissoc: " + str(fails))
print("\nMCH basis")
print("     State 1: " + str(s21) + "       State 2: " + str(s22) + "       No dissoc: " + str(fails))
print()
