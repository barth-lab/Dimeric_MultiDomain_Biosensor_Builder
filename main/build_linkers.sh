#!/bin/bash

# https://stackoverflow.com/questions/14447406/bash-shell-script-check-for-a-flag-and-grab-its-value
# https://archive.is/TRzn4
h="true"
r=""
o=0
l=0
while getopts ":R:C:h:d:r:o:l:" opt; do
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
    d)
      set -f # disable glob
      IFS=' ' # split on space characters
      d=($OPTARG) ;; # dimerisation sites
    r)
      r=$OPTARG # what round are we building from (e.g. if we needed intermediate stages to include domains)
      # this needs to be _2 or whatever (not just 2)
      ;;
    o)
      set -f
      IFS=' '
      o=($OPTARG) # reorder the domains into what order based on previous input (basically this needs to address missing domains, creating two unique chains)
      # based on all_verbose.fasta in the last phase, what is the needed order here? This will be larger than those number of domains with dimerisation domains
      ;; 
    l)
      l=$OPTARG # Position of the ligand in terms of the original all_verbose.fasta (most likely 1)
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

# normal usage
# dimer domains needs to now be just the domain that is still treated as a dimer (i.e. two repeating sequneces)
# bash /data/domain_construction/domain_assembly_constraints/main/build_linkers.sh -R /data/rosetta20_glis -C 7.29 -d "3" -r _2 -o "4 5 6 3 1 2"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# prep files
I=input_loop
mkdir $I
O=output_loop
mkdir $O

# permit wildcard use
set +o noglob

# grab all cluster centers and put them in the input loop folder
cp ./output_scaffold${r}/${C}/c.*.0.pdb ./${I}
cp ./input_scaffold${r}/all_verbose.fasta ./${I} # - use this as a basis to reorder PDB
cp ./input_scaffold${r}/cst ./${I} # cst contain resid positions where we need to rebuild between
cp ./input_scaffold${r}/frags* ./${I}

# Add in residues at the correct positions for rebuilding
cnt=-1
for i in ./${I}/c.*.0.pdb; do
    python ${SCRIPT_DIR}/get_resid_reorder.py -s ${i} -d "${d[@]}"  -f ./input_scaffold${r}/all_verbose.fasta -o "${o[@]}" -l $l
    tac ${i} | awk '/TER/ {if (f) next; f=1}1' | tac > tmp # remove duplicate TER lines
    mv tmp ${i}
    cnt=$(expr $cnt + 1)
done
# get_resid_reorder.py also creates the loopfile needed
mv loopfile ./${I}
mv new_verbose.fasta ./${I}

# remove unncessary TER lines from insertion as Rosetta doesn't like this

if [ "$h" = "true" ] ; then
    echo 'Now run mp_buildloops.slurm HPC command in hpc_run/ folder of main/'
    echo "Change n in the mp_loopbuild file to $cnt"
else
    # note this is meant more as a check, and will only run loop model on the lowest energy structure in input
    $R/source/bin/loopmodel.linuxgccdebug \
        -database $R/database/  \
        -in:file:s $I/c.0.0.pdb \
        -loop_file $I/loopfile \
        -loops:frag_sizes 9 3 \
        -loops:frag_files input_scaffold${r}/frags.200.9mers input_scaffold${r}/frags.200.3mers \
        -loops:remodel quick_ccd \
        -loops:refine refine_ccd \
        -loops:relax no \
        -constant_seed \
        -nstruct 1 \
        -out:suffix _loop \
        -out:path:all $O \
        -out:pdb \
	-ignore_zero_occupancy false \
        > $O/log \
        2> $O/err
fi

