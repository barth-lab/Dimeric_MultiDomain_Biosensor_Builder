#!/bin/bash

# https://stackoverflow.com/questions/14447406/bash-shell-script-check-for-a-flag-and-grab-its-value
# https://archive.is/TRzn4
while getopts ":R:T:d:" opt; do
  case $opt in
    R)
      #echo "-R was triggered, Parameter: $OPTARG" >&2
      R=$OPTARG # rosetta location
      ;;
    T) 
      TM=$OPTARG # which input domain corresponds to the TM region
      ;;
    s)
      # allow importing of multiple domains https://unix.stackexchange.com/questions/164259/provide-two-arguments-to-one-option-using-getopts 
      set -f # disable glob
      IFS=' ' # split on space characters
      domains=($OPTARG) ;; # use the split+glob operator
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

O=output_scaffold
mkdir $O

# create fragment files based on fasta
# https://www.rosettacommons.org/docs/latest/application_documentation/utilities/app-fragment-picker
$R/source/bin/fragment_picker.static.linuxgccrelease \
    –database $R/main/database/  \
    -in:file:vall
    -in:file:fasta
    frags:scoring:config
    frags:frag_sizes
    frags:picking:query_pos # linker positions from fasta file!

# now run domain assembly
$R/source/bin/mp_domain_assembly.linuxgccrelease \
    -database $R/main/database/  \
    -in:file:fasta input_scaffold/all.fasta \
    -in:file:frag3 input_scaffold/03.frags \
    -in:file:frag9 input_scaffols/09.frags \
    -mp:assembly:TM_pose_number ${TM} \
    -mp:assembly:poses $(for i in "${domains[@]}"; do echo ${i}; done) \
    -rebuild_disulf false \
    -detect_disulf false \
    -constraints::cst_fa_file input_scaffold/cst \
    -constraints:cst_fa_weight 1 \
    -constant_seed \
    -nstruct 1 \
    -out:path:all $O \
    -out:pdb \
    > $O/log \
    2> $O/err
