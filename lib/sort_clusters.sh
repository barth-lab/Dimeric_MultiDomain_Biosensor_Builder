#!/bin/bash

round() {
  printf "%.${2}f" "${1}"
}

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

size_var=(0.5 0.533 0.567 0.6 0.633 0.667 0.7 0.733 0.767 0.8 0.833 0.867)

big_radius=${1}

loc=$(pwd)
for s in ${size_var[@]}; do
    rad=$(bc <<< "${big_radius} * ${s}")
    num=$(round $rad 2)
    cd ${num}

    # extract cluster data between certain values
    #sed -n '/protocols.cluster: ---------- Summary ---------------------------------/{:a;n;/protocols.cluster: ----------------------------------------------------/b;p;ba}' clout > cluster_summary.dat

    # now check the stats in each cluster
    python ${SCRIPT_DIR}/sort_clusters.py 
    cd ${loc}
done
