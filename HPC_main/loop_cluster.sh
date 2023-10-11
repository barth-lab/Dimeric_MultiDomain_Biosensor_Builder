#!/bin/bash

round() {
  printf "%.${2}f" "${1}"
}
size_var=(0.5 0.533 0.567 0.6 0.633 0.667 0.7 0.733 0.767 0.8 0.833 0.867)

in_file=$1
big_radius=$2

loc=$(pwd)

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

for s in ${size_var[@]}; do
    rad=$(bc <<< "${big_radius} * ${s}")
    num=$(round $rad 2)
    mkdir ${num}
    cd ${num}
    sbatch ${SCRIPT_DIR}/mp_cluster.slurm ${loc}/${in_file} ${num}
    cd ${loc}
done


