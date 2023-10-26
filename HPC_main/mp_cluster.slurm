#!/bin/bash
#SBATCH --time 12:00:00 -n 1
#SBATCH --partition=serial
#SBATCH --job-name=cluster

# run on fidis

R=/location/to/rosetta
in_file=$1
radius=$2

echo ${current}

$R/source/bin/cluster.linuxgccrelease \
        -database $R/database/  \
        -silent_read_through_errors \
        -in:file:fullatom \
        -in:file:silent $1 \
        -cluster:radius $radius \
        -run:shuffle \
        -cluster:population_weight 1 \
        -cluster:sort_groups_by_energy \
        -cluster:write_centers \
        > clout \
        2> clerr

#-out:path:all output_scaffold/$radius/ \
