"""
Rosetta strangely doesn't like biobox written PDB file input - so sort here
"""

import sys
filename=sys.argv[1]

with open(filename, "r") as f:
    lines = f.readlines()

print(lines)

new_pdb=sys.argv[1]
with open(new_pdb, "w") as f:
    for line in lines:
        if line[:3] == "TER" or line[:3] == "END" or line[:5] == "MODEL":
            pass
        else:
            f.write(line)

