"""
Check if final output structures from loop building, if any, are decent quality
"""

import biobox as bb
import numpy as np
import sys, os

def check_structure(M):
    """
    Check coordinates of M
    """
    pts = M.points
    if np.any(pts > 900):
        return False
    else:
        return True

directory=str(sys.argv[1])

for file in os.listdir(directory):
     filename = os.fsdecode(file)
     if filename.endswith(".pdb"):
         M = bb.Molecule(directory + "/" + filename)
         check = check_structure(M)
         if check:
             print(filename)
