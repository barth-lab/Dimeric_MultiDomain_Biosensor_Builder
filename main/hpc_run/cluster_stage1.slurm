#!/bin/bash
#SBATCH --time 10:00:00 -n 1
#SBATCH --mem 48G
#SBATCH --partition=serial
#SBATCH --job-name=cluster

# run on fidis

R=/home/dclw/rosetta20_glis/
in_file=../filtered.silent
radius=$1

$R/source/bin/cluster.linuxgccdebug \
        -database $R/database/  \
        -silent_read_through_errors \
        -in:file:fullatom \
        -in:file:silent $in_file \
        -cluster:radius $radius \
        -run:shuffle \
        -cluster:population_weight 1 \
        -cluster:sort_groups_by_energy \
        -cluster:write_centers \
        > clout \
        2> clerr

#-out:path:all output_scaffold/$radius/ \
