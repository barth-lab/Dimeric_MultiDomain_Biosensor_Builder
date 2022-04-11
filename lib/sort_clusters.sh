#!/bin/bash

round() {
  printf "%.${2}f" "${1}"
}

size_var=(0.5 0.533 0.567 0.6 0.633 0.667 0.7 0.733 0.767 0.8)

big_radius=11.825

loc=$(pwd)
#for i in 3.84 4.14 4.44 4.74 5.64 6.01 6.39 6.76 7.14 7.52 7.89 8.26 8.65 9.02 9.39 9.77; do
for s in ${size_var[@]}; do
    rad=$(bc <<< "${big_radius} * ${s}")
    num=$(round $rad 2)
    echo ${num}
    cd ${num}

    # extract cluster data between certain values
    #sed -n '/protocols.cluster: ---------- Summary ---------------------------------/{:a;n;/protocols.cluster: ----------------------------------------------------/b;p;ba}' clout > cluster_summary.dat

    # now check the stats in each cluster
    # I then count how many clusters I get and how many structures there are in the biggest cluster for each radius, and chose the one which seems the most reasonable, i.e. at least 20 clusters, between 10% and 20% of all structures in the largest cluster, at least 1-5% in the smallest cluster
    python /data/domain_construction/domain_assembly_constraints/lib/sort_clusters.py #cluster_summary.dat
    cd ${loc}
done
