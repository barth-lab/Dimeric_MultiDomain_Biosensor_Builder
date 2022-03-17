#!/bin/bash

for i in 3.30 3.52 3.74 3.90 4.18 4.40 4.60 4.84 5.06 5.20; do
    cd ${i}
    sbatch /scratch/dclw/domain_assembly/domain_assembly_constraints/main/hpc_run/cluster_stage1.slurm ${i}
    cd ..
done
