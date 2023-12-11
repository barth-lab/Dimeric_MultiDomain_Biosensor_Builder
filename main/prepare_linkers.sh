#!/bin/bash

S="2"
o=0
l=0
while getopts ":C:d:S:o:l:" opt; do
  case $opt in
    C)
      # cluster folder to select (extract all centers from this folder)
      C=$OPTARG
      ;;
    d)
      set -f # disable glob
      IFS=' ' # split on space characters
      d=($OPTARG) ;; # dimerisation sites in terms of input all_verbose.fasta
    S)
      r=$OPTARG # what round are we building from (e.g. if we needed intermediate stages to include domains)
      # this needs to be 2 or whatever
      ;;
    o)
      set -f
      IFS=' '
      o=($OPTARG) # reorder the domains into what order based on previous input (basically this needs to address missing domains, creating two unique chains)
      # based on all_verbose.fasta in the last phase, what is the needed order here? This will be larger than those number of domains with dimerisation domains
      ;; 
    l)
      l=$OPTARG # Position of the ligand in terms of the previous all_verbose.fasta (most likely 1). Only include if there is a ligand of course!
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
# bash /data/domain_construction/domain_assembly_constraints/main/build_linkers.sh -R /data/rosetta20_glis -C 6.25 -d "1 2 3" -o "1 4 5 2 3 6" -l 1
# SCF extra linker final stage:
# bash /data/domain_construction/domain_assembly_constraints/main/build_linkers.sh -R /data/rosetta20_glis -C 8.09 -d "3 4" -o "3 6 7 8 4 5 1 2" -l 3 -r _2

# This script needs to be run from the run folder

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# prep files
I=input_loop
mkdir $I
O=output_loop
mkdir $O

# permit wildcard use
set +o noglob

# grab all cluster centers and put them in the input loop folder
cp ./output_scaffold_${S}/${C}/c.*.0.pdb ./${I}
cp ./input_scaffold_${S}/all_verbose.fasta ./${I} # - use this as a basis to reorder PDB
cp ./input_scaffold_${S}/cst ./${I} # cst contain resid positions where we need to rebuild between
cp ./input_scaffold_${S}/frags* ./${I}

# Add in residues at the correct positions for rebuilding
for i in ./${I}/c.*.0.pdb; do
    python ${SCRIPT_DIR}/get_resid_reorder.py -s ${i} -d "${d[@]}"  -f ./input_scaffold${r}/all_verbose.fasta -o "${o[@]}" -l $l
    tac ${i} | awk '/TER/ {if (f) next; f=1}1' | tac > tmp # remove duplicate TER lines
    mv tmp ${i}
done
# get_resid_reorder.py also creates the loopfile needed
mv loopfile ./${I}
mv new_verbose.fasta ./${I}
