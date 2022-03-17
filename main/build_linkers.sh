#!/bin/bash

# https://stackoverflow.com/questions/14447406/bash-shell-script-check-for-a-flag-and-grab-its-value
# https://archive.is/TRzn4
h=true
while getopts ":R:C:h:" opt; do
  case $opt in
    R)
      #echo "-R was triggered, Parameter: $OPTARG" >&2
      R=$OPTARG # rosetta location
      ;;
    C)
      # cluster folder to select (extract all centers from this folder)
      C=$OPTARG
      ;;
    h)
      # run as HPC (i.e. don't run loopmodelling here, just prep)
      # default is 1 (i.e. don't run loopmodel here)
      h=$OPTARG
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

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# prep files
I=input_loop
mkdir $I
O=output_loop
mkdir $O

# grab all cluster centers and put them in the input loop folder
cp ./output_scaffold/${C}/c.*.0.pdb ${I}
cp ./input_scaffold/all_verbose.fasta ${I} # - use this as a basis to reorder PDB
cp ./input_scaffold/cst ${I} # cst contain resid positions where we need to rebuild between

# Add in residues at the correct positions for rebuilding

# now create the loopfile


if [ "$h" = true ] ; then
    echo 'Now run loopmodel HPC command in hpc_run/ folder of main/'
else
    # note this is meant more as a check, and will only run loop model on the lowest energy structure in input
    $R/source/bin/loopmodel.linuxgccdebug \
        -database $R/database/  \
        -in:file:s $I/c.0.0.pdb \
        -loop_file $I/loopfile \
        -loops:frag_sizes 9 3 \
        -loops:frag_files input_scaffold/frags.200.9mers input_scaffold/frags.200.3mers \
        -loops:remodel quick_ccd \
        -loops:refine refine_ccd \
        -loops:relax no \
        -constant_seed \
        -nstruct 1 \
        -out:suffix _loop \
        -out:path:all $O \
        -out:pdb \
        > $O/log \
        2> $O/er
fi

