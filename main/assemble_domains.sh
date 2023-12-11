#!/bin/bash
# Construction always needs to occur down (i.e. towards cytoplasm)

S=2
# define flags
while getopts ":S:s:d:C:p:" opt; do
  case $opt in
    S)
      S=$OPTARG # which stage of assembly are we on (default is 2)
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
      p=$OPTARG # current (i.e. prior to assembly) position of the domain(s) being added - upgrade to include more domains
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

# make run folder - remove for testing to avoid repeat calculation
# stage S of assembly
touch input_scaffold_${S}
rm -r input_scaffold_${S}
mkdir input_scaffold_${S}
cd input_scaffold_${S}/

# if we are on stage S, what was the previous stage?
if [ $S = 2 ]; then
    stage0=""
else
    ((num= S - 1))
    stage0="_${num}"
fi

# need to change in case two domains are added
name=$(echo "${domains[0]}" | cut -d'.' -f1)

set +o noglob
cp ../input_scaffold${stage0}/${name}_cut.pdb .
cp ../output_scaffold${stage0}/${C}/c.*.0.pdb .
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
