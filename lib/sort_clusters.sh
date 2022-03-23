#!/bin/bash

loc=$(pwd)
#for i in 3.84 4.14 4.44 4.74 5.64 6.01 6.39 6.76 7.14 7.52 7.89 8.26 8.65 9.02 9.39 9.77; do
for i in 3.34 3.64 3.94 4.24 4.54 4.84 5.14 5.45 5.74 6.05 6.35 6.65 6.96 7.26 7.56 7.87; do
    echo ${i}
    cd ${i}

    # extract cluster data between certain values
    #sed -n '/protocols.cluster: ---------- Summary ---------------------------------/{:a;n;/protocols.cluster: ----------------------------------------------------/b;p;ba}' clout > cluster_summary.dat

    # now check the stats in each cluster
    # I then count how many clusters I get and how many structures there are in the biggest cluster for each radius, and chose the one which seems the most reasonable, i.e. at least 20 clusters, between 10% and 20% of all structures in the largest cluster, at least 1-5% in the smallest cluster
    python /data/domain_construction/domain_assembly_constraints/lib/sort_clusters.py #cluster_summary.dat
    cd ${loc}
done
