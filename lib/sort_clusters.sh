#!/bin/bash

# extract cluster data between certain values
sed -n '/protocols.cluster: ---------- Summary ---------------------------------/{:a;n;/protocols.cluster: ----------------------------------------------------/b;p;ba}' clout > cluster_summary.dat

# now check the stats in each cluster
# I then count how many clusters I get and how many structures there are in the biggest cluster for each radius, and chose the one which seems the most reasonable, i.e. at least 20 clusters, between 10% and 20% of all structures in the largest cluster, at least 1-5% in the smallest cluster
python /data/domain_construction/domain_assembly_constraints/lib/sort_clusters.py cluster_summary.dat
