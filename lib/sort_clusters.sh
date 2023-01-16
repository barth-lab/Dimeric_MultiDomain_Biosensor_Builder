#!/bin/bash

round() {
  printf "%.${2}f" "${1}"
}

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

size_var=(0.5 0.533 0.567 0.6 0.633 0.667 0.7 0.733 0.767 0.8)

big_radius=${1}

loc=$(pwd)
for s in ${size_var[@]}; do
    rad=$(bc <<< "${big_radius} * ${s}")
    num=$(round $rad 2)
    echo ${num}
    cd ${num}

    # extract cluster data between certain values
    #sed -n '/protocols.cluster: ---------- Summary ---------------------------------/{:a;n;/protocols.cluster: ----------------------------------------------------/b;p;ba}' clout > cluster_summary.dat

    # now check the stats in each cluster
    # I then count how many clusters I get and how many structures there are in the biggest cluster for each radius, and chose the one which seems the most reasonable, i.e. at least 20 clusters, between 10% and 20% of all structures in the largest cluster, at least 1-5% in the smallest cluster
    python ${SCRIPT_DIR}/sort_clusters.py #cluster_summary.dat
    cd ${loc}
done
