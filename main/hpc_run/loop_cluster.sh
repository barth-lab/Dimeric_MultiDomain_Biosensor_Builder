#!/bin/bash

round() {
  printf "%.${2}f" "${1}"
}
size_var=(0.5 0.533 0.567 0.6 0.633 0.667 0.7 0.733 0.767 0.8)

big_radius=14.575

loc=$(pwd)

for s in ${size_var[@]}; do
    rad=$(bc <<< "${big_radius} * ${s}")
    num=$(round $rad 2)
    mkdir ${num}
    cd ${num}
    sbatch /scratch/dclw/domain_assembly/domain_assembly_constraints/main/hpc_run/mp_cluster.slurm ${num}
    cd ${loc}
done


