#!/bin/bash

# $1 rosetta loc
# $2 silent files
# $3 radius (initial -1 for guessing)
# $4 location of output files

$1/source/bin/cluster.linuxgccdebug \
        -database $1/database/  \
        -silent_read_through_errors \
        -in:file:fullatom \
        -in:file:silent $2 \
        -cluster:radius $3 \
        -run:shuffle \
        -cluster:population_weight 1 \
        -cluster:sort_groups_by_energy \
        -cluster:write_centers \
        > $4/clout \
        2> $4/clerr
