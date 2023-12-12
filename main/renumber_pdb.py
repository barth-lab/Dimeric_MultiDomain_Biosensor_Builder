import biobox as bb
import sys, os
import numpy as np

A = bb.Molecule(sys.argv[1])

# renumber from 1 to x, all resid. Reset chain ID to A also
CA_idx = np.asarray(A.atomselect("*", "*", "CA", get_index=True)[1])
resnum = np.asarray(A.data['resid'][CA_idx])
# chain for each resid
chains = np.asarray(A.data['chain'][CA_idx])


# start residue numbering from 1. Change when chain break occurs (in file, not in structure)
res_count = 1
for cnt, val in enumerate(CA_idx):
    # maximum AA length is 27 (tryp with hydrogens), set greater than 30 as threashold
    # full residue index set
    full_res = A.atomselect(chains[cnt], [resnum[cnt]], "*", get_index=True)[1]

    # now remove residues that have similar properties, but are not the same
    full_res = np.asarray([x for x in full_res - val if np.abs(x) <= 30]) + val

    # now renumber
    A.data.loc[full_res, "resid"] = res_count

    res_count += 1

A.data["chain"] = "A"

A.write_pdb(sys.argv[1])
