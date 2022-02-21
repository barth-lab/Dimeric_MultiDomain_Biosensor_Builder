#!/bin/bash

# https://stackoverflow.com/questions/14447406/bash-shell-script-check-for-a-flag-and-grab-its-value
# https://archive.is/TRzn4
while getopts ":R:T:s:" opt; do
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

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

O=output_scaffold
mkdir $O

# don't run frag picker unless necessary
if [ ! -f input_scaffold/frags.200.3mers ] || [ ! -f input_scaffold/frags.200.9mers ]; then 
    # get linker positions for where fragment picker needs to occur (include neighbouring sequences)
    linker_pos=$(python $SCRIPT_DIR/frag_positions.py)

    # create fragment files based on fasta
    # https://www.rosettacommons.org/docs/latest/application_documentation/utilities/app-fragment-picker
    # eventually this needs to be improved because it has no secondary structure stuff in it (PSIPRED needed!!)
    $R/main/source/bin/fragment_picker.static.linuxgccrelease \
        -database $R/main/database/  \
        -in:file:vall $R/main/tools/fragment_tools/vall.jul19.2011.gz \
        -in:file:fasta input_scaffold/all.fasta \
        -frags:scoring:config $SCRIPT_DIR/../lib/simple.wghts \
        -frags:frag_sizes 3 9 \
        -frags:picking:query_pos $(echo ${linker_pos}) # linker positions from fasta file!

    mv frags.200.3mers frags.200.9mers input_scaffold/    
fi

cd input_scaffold/
set +o noglob
for f in *.pdb; do
    python $SCRIPT_DIR/prepare_input_pdb.py ${f}
done
cd ../

# now run domain assembly
# need to resolve static etc. later
$R/main/source/bin/mp_domain_assembly.static.linuxgccrelease \
    -database $R/main/database/  \
    -in:file:fasta input_scaffold/all.fasta \
    -in:file:frag3 input_scaffold/frags.200.3mers \
    -in:file:frag9 input_scaffold/frags.200.9mers \
    -mp:assembly:TM_pose_number ${TM} \
    -mp:assembly:poses "${domains[@]}" \
    -rebuild_disulf false \
    -detect_disulf false \
    -constraints::cst_fa_file input_scaffold/cst \
    -constraints:cst_fa_weight 1 \
    -constant_seed \
    -nstruct 1 \
    -out:path:all $O \
    -out:pdb \
    -ignore_zero_occupancy false \
    > $O/log \
    2> $O/err
