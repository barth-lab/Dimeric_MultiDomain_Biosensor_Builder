"""
Check stats in import cluster summary file 
"""

import numpy as np
import sys, os
import glob

# criteria:
# I then count how many clusters I get and how many structures there are in the biggest cluster for each radius, and chose the one which seems the most reasonable, i.e. at least 20 clusters, between 10% and 20% of all structures in the largest cluster, at least 1-5% in the smallest cluster

def check_stats(directory):
    """
    Check the clustering statistics in filename
    """
    C = [] # each cluster is appended to this global array
    # go through files in directoy, assuming naming convention of c.C.N.pdb (C = cluster no., N = number in cluster
    # in increasing energy)
    counter = True
    cnt = 0
    total_len = 0
    while counter: # don't know how many cluster files there are
        
        ClusterCounter = len(glob.glob1(directory,"c.%i.*"%(cnt)))
        if ClusterCounter == 0:
            counter = False
        else:
            cnt += 1
            C.append(ClusterCounter)
            total_len += ClusterCounter

    #with open(filename, "r") as f:
    #    lines = f.readlines()
    #    c_len = 0 # current cluster length
    #    total_len = 0
    #    for line in lines:
    #       if line[19:26] == "Cluster":
    #           # reset count
    #           C.append(c_len)
    #           c_len = 0
    #       else:
    #           c_len += 1
    #           total_len += 1

    #C = C[1:] # first appended length is a dummy (as first line in file contains Cluster)
    ##### Now print stats ####
    print("#######")
    print("")
    print(">> Number of clusters: %i"%(len(C)))
    print(">> Number of structures in largest cluster: %.2f %%"%(100 * float(C[0]) / total_len))
    print(">> Number of structures in smallest cluster: %.2f %%"%(100 *float(C[-1]) / total_len))
    print("")
    print("#######")

# c.0.0.pdb is the lowest energy model from the first cluster, etc.
#filename = sys.argv[1]
path=os.getcwd()
check_stats(path)
