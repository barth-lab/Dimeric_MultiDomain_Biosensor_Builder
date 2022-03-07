#!/bin/bash

# https://stackoverflow.com/questions/14447406/bash-shell-script-check-for-a-flag-and-grab-its-value
# https://archive.is/TRzn4
while getopts ":R:" opt; do
  case $opt in
    R)
      #echo "-R was triggered, Parameter: $OPTARG" >&2
      R=$OPTARG # rosetta location
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

echo $R

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# don't run frag picker unless necessary
if [ ! -f input_scaffold/frags.200.3mers ] || [ ! -f input_scaffold/frags.200.9mers ]; then 
    # get linker positions for where fragment picker needs to occur (include neighbouring sequences)
    linker_pos=$(python $SCRIPT_DIR/../frag_positions.py)

    # create fragment files based on fasta
    # https://www.rosettacommons.org/docs/latest/application_documentation/utilities/app-fragment-picker
    # eventually this needs to be improved because it has no secondary structure stuff in it (PSIPRED needed!!)
    $R/source/bin/fragment_picker.linuxgccdebug \
        -database $R/database/  \
        -in:file:vall $R/tools/fragment_tools/vall.jul19.2011.gz \
        -in:file:fasta input_scaffold/all.fasta \
        -frags:scoring:config $SCRIPT_DIR/../../lib/simple.wghts \
        -frags:frag_sizes 3 9 \
        -frags:picking:query_pos $(echo ${linker_pos}) # linker positions from fasta file!

    mv frags.200.3mers frags.200.9mers input_scaffold/    
fi

cd input_scaffold/
set +o noglob
for f in *.pdb; do
    python $SCRIPT_DIR/../prepare_input_pdb.py ${f}
done
cd ../
