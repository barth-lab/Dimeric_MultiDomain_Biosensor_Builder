"""
Simple script to split the chains of a molecule such that idealise/relax can run correctly
"""
import biobox as bb
import sys, os
import numpy as np

loop=True

if loop:
    # because of loop insertions you can have weird chain breaks, so use different method to assign chains
    M = bb.Molecule(sys.argv[1])
    no_res = np.sum(M.data["name"] == "CA")

    # currently resid are 1 to no_res from output - so take halves
    resid_A = np.arange(1, no_res/2+1).astype(int)
    resid_B = np.arange(no_res/2+1, no_res+1).astype(int)
    # +1 because of resid numbering

    idx_A = M.atomselect("*", resid_A, "*", get_index=True)[1]
    idx_B = M.atomselect("*", resid_B, "*", get_index=True)[1]

    #M.data["chain"].iloc[idx_B] = "B"
    M.data["chain"].iloc[idx_B] = "B"

    M.write_pdb(sys.argv[1])
else:
    M = bb.Molecule(sys.argv[1])
    M.guess_chain_split()
    M.write_pdb(sys.argv[1])
