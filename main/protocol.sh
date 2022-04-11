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
h="true"

# define flags
while getopts ":R:s:T:l:d:x:a:h:" opt; do
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
    x)
      x="-x $OPTARG" # extra_linker file location
      ;;
    a)
      set -f
      IFS=' '
      a=($OPTARG) 
      a="-a "${a[@]}"";; # avoid linker cutting domains
    h)
      h=$OPTARG # is the plan to run this on a HPC?
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

# current script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# make run folder - remove for testing to avoid repeat calculations, TODO: convert this to checkpoint so we don't repeat from start each time
D=`date +"%s"`
mkdir run${D}
cp "${domains[@]}" ${l} ${x} run${D}/
cd run${D}
mkdir input_scaffold
#cd run1645438828/

# for each input domain, check if it's a dimer (to see where crossing points need to be for design) - note, this need to avoid checking for ligand!
# FOr now this can be done by seeing if its a homodimer (how to do heterodimers later?)
# Then we need to remove the linkers (but include within fasta)
# then we need to build the fragment database
# THen we run domain assembly

# OTHER NEEDED INPUT FILES:
# 1) fragment files (for domain assembly and linker modelling) - can be built with fasta file created by prepare_scaffold
# 2) fasta file of domains and, eventually, chimera
# 3) Span file for the TM domain

# two protocols - one to run on cluster, one without?

# read input domains from a file is probably easiest
# Starting scaffold needs the start and end points included (i.e. EC head to CT region) to build inial scaffold. This is based on the input domains
${SCRIPT_DIR}/prepare_scaffold.py -s "${domains[@]}" -d "${d[@]}" -l ${l} ${a} ${x}

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
echo $domains_str

if [ "$h" == "true" ] ; then
    echo "Now use the fasta file to create your own fasta fragments (e.g. using Robertta)"
    echo "before running the mp_assembly on a HPC"    
else
    # run round of first stage of assembly
    bash ${SCRIPT_DIR}/mp_assembly_stage1.sh -R $R -T ${TM} -s "$(echo $domains_str)"
fi

# need to cut the linker regions during assembly as they'll be readded later - however we don't want to apply this to D3 as this is unstructured already

# run domain assembly with constraints - building from top to bottom (EC to CT)

# cluster

# run linker design

# cluster?

# select from miminum energy, and repear
