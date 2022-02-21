"""
Get linker positions for fragment picker (including neighbouring resid for smoother fragments)
"""

import numpy as np
import sys, os

with open("input_scaffold/all_verbose.fasta", "r") as f:
#with open("all_verbose.fasta", "r") as f:
    lines = f.readlines()

linker_pos = []
fasta_cnt = 1 # 1 because of rosetta indexing, would be 0 if python
prepare_linker = False
for line in lines:
    if line[:7] == ">linker": # linker line incoming
        prepare_linker = True
    elif line[0] == ">": # no linker line incoming
        prepare_linker = False
    elif prepare_linker:
        linker_pos.append(np.arange(fasta_cnt-1, fasta_cnt + len(line))) # get peripheral areas also!
        fasta_cnt += len(line.strip()) # in case there is whitespace, strip
    else:
        fasta_cnt += len(line.strip())

linker_pos = np.concatenate(linker_pos).ravel()

print(' '.join(linker_pos.astype(str)))

### testing only ###
#with open("all.fasta", "r") as f:
#    fasta = f.readlines()[1]

#for l in linker_pos:
#    print(fasta[l])