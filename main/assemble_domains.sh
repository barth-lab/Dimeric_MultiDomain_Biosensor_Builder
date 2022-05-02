#!/bin/bash

# author: Lucas Rudden

# Here I want to write out the entire protocol of the modified approach to domain assembly
# Assembly procotol is essentially domain assemble (with constraints), clustered and then taken further for assembly rounds
# Construction always needs to occur down (i.e. towards cytoplasm)

# Note that static vs mpi vs nothing in linux rosetta commands

# https://stackoverflow.com/questions/14447406/bash-shell-script-check-for-a-flag-and-grab-its-value
# https://archive.is/TRzn4

# optional flags need default values to avoid errors
h="true"

# define flags
while getopts ":R:s:T:d:C:p:" opt; do
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
    d)
      set -f # disable glob
      IFS=' ' # split on space characters
      d=($OPTARG) ;; # dimerisation sites
    C) 
      C=$OPTARG # cluster center to continue assembly from
      ;;
    p) 
      p=$OPTARG # position of the domain(s) being added - upgrade to include more domains
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

# example of this script being run
# bash /data/domain_construction/domain_assembly_constraints/main/assemble_domains.sh -R /data/rosetta20_glis -T 3 -d "1 3" -p 2 -C 7.89 -s "D2.pdb c.0.0.pdb"

# current script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# make run folder - remove for testing to avoid repeat calculations, TODO: convert this to checkpoint so we don't repeat from start each time
# stage 2 of assembly
mkdir input_scaffold_2
cd input_scaffold_2/

# need to change in case two domains are added
name=$(echo "${domains[0]}" | cut -d'.' -f1)

set +o noglob
cp ../input_scaffold/${name}_cut.pdb .
cp ../output_scaffold/${C}/c.*.0.pdb .
# copy frag files for linker insertion during domain assembly
cp ../input_scaffold/*frags* .
# any of the cluster centers will work for the constraint and fasta files, since sequence wise they are identical
# e.g. python /data/domain_construction/domain_assembly_constraints/main/assemble_domains.py -s D2.pdb c.0.0.pdb -d 1 3 --add_domains 2

# for each cluster center, we need to restructure the metadata
for s in c.*.0.pdb; do
    ${SCRIPT_DIR}/assemble_domains.py -s ${name}_cut.pdb ${s} -d "${d[@]}" --add_domains ${p}
    tac ${s} | awk '/TER/ {if (f) next; f=1}1' | tac > tmp # remove duplicate TER lines
    mv tmp ${s}
done

if [ "$h" == "true" ] ; then
    echo "now run mp_assembly_stage2.sh on the HPC with parameters"
else
    # run round of first stage of assembly
    bash ${SCRIPT_DIR}/mp_assembly_stage2.sh -R $R -T ${TM} -s "$(echo $domains_str)"
fi

# need to cut the linker regions during assembly as they'll be readded later - however we don't want to apply this to D3 as this is unstructured already
# run domain assembly with constraints - building from top to bottom (EC to CT)
# cluster
# run linker design
# cluster
# select from miminum energy, and repear
