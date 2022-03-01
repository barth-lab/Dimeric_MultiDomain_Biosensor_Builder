#!/bin/bash
#SBATCH --time 05:00:00 -n 1
#SBATCH --mem 48G
#SBATCH --partition=serial

# run on fidis

R=/home/dclw/rosetta20_glis/
in_file=
radius=

$R/source/bin/cluster.linuxgccdebug \
        -database $R/database/  \
        -silent_read_through_errors \
        -in:file:fullatom \
        -in:file:silent $1 \
        -cluster:radius $2 \
        -run:shuffle \
        -cluster:population_weight 1 \
        -cluster:sort_groups_by_energy \
        -cluster:write_centers \
        > clout \
        2> clerr
