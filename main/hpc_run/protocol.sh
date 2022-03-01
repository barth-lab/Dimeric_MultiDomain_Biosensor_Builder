#!/bin/bash

# author: Lucas Rudden

# Here I want to write out the entire protocol of the modified approach to domain assembly
# Assembly procotol is essentially domain assemble (with constraints), clustered and then taken further for assembly rounds
# Construction always needs to occur down (i.e. towards cytoplasm)

# Note that static vs mpi vs nothing in linux rosetta commands

# https://stackoverflow.com/questions/14447406/bash-shell-script-check-for-a-flag-and-grab-its-value
# https://archive.is/TRzn4

# optional flags need default values to avoid errors
x=""
a=""

# define flags
while getopts ":R:s:T:l:d:x:a:S:" opt; do
  case $opt in
    R)
      R=$OPTARG # rosetta location
      ;;
    T)
      TM=$OPTARG # which input domain corresponds to the TM region
      ;;
    s)
      # allow importing of multiple domains https://unix.stackexchange.com/questions/164259/provide-two-arguments-to-one-option-using-getopts 
      set -f # disable glob
      IFS=' ' # split on space characters
      domains=($OPTARG) ;; # domain pdbs
    l) 
      l=$OPTARG # linker file (location of where to cut linkers)
      ;;
    d)
      set -f # disable glob
      IFS=' ' # split on space characters
      d=($OPTARG) ;; # dimerisation sites
    S) 
      # stage of the procotol to run (e.g. 1, 2...) split so we can run on HPC easier
      S=$OPTARG
      ;;
    x)
      x="-x $OPTARG" # extra_linker file location
      ;;
    a)
      set -f
      IFS=' '
      a=($OPTARG) 
      a="-a "${a[@]}"";; # avoid linker cutting domains
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

# current script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

if [ $S == 1 ]; then

    # make run folder - remove for testing to avoid repeat calculations, TODO: convert this to checkpoint so we don't repeat from start each time
    D=`date +"%s"`
    mkdir run${D}
    cp "${domains[@]}" ${l} ${x} run${D}/
    cd run${D}
    mkdir input_scaffold

    # read input domains from a file is probably easiest
    # Starting scaffold needs the start and end points included (i.e. EC head to CT region) to build inial scaffold. This is based on the input domains
    ${SCRIPT_DIR}/../prepare_scaffold.py -s "${domains[@]}" -d "${d[@]}" -l ${l} ${a} ${x}

    mv cst input_scaffold/ 
    find -name "*_cut.pdb" -type f -exec mv -t input_scaffold/ {} +
    find -name "*fasta" -type f -exec mv -t input_scaffold/ {} +

    # rename domains to respect the addition of the _cut suffix
    domains_cut=()
    for i in "${domains[@]}"; do
        name=$(echo $i | cut -d'.' -f1)
        domains_cut+=("input_scaffold/${name}_cut.pdb")
    done

    printf -v domains_str ' %s' "${domains_cut[@]}"

    ${SCRIPT_DIR}/create_frags.sh $R

    echo "########  Protocol stage 1 complete #########"
    echo ""
    echo ">> Check input_scaffold - fragments created, fastas, and constraints prepared."
    echo ">> Now run on the HPC mp_assembly_stage1.slurm with the following variables: "
    echo "TM=$T"
    echo "R=/home/dclw/rosetta20_glis/" # rosetta location with assembly protocol
    echo "domains=$domains_str"
    echo ""
    echo "##############################################"
fi



# run domain assembly with constraints - building from top to bottom (EC to CT)

# cluster

# run linker design

# cluster?

# select from miminum energy, and repear
