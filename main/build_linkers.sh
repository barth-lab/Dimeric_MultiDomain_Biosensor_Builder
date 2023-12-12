#!/bin/bash

N=1
# define flags
while getopts ":R:N:" opt; do
  case $opt in
    R)
      R=$OPTARG # Location of Rosetta
      ;;
    N)
      N=$OPTARG # number of desired output models
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

I=input_loop
O=output_loop

# loop through all cluster centres
for i in $I/c*.pdb; do
    $R/source/bin/loopmodel.linuxgccrelease \
        -database $R/database/  \
        -in:file:s ${i} \
        -loop_file $I/loopfile \
        -loops:frag_sizes 9 3 \
        -loops:frag_files input_scaffold/frags.200.9mers input_scaffold/frags.200.3mers \
        -loops:remodel quick_ccd \
        -loops:refine refine_ccd \
        -loops:relax no \
        -constant_seed \
        -nstruct $N \
        -out:suffix _loop \
        -out:path:all $O \
        -out:pdb \
        -ignore_zero_occupancy false \
        > $O/log \
        2> $O/err
done

