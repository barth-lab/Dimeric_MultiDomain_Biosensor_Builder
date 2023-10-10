#!/bin/bash
# author: Lucas Rudden

# prepare the initial scaffold for running the domain assembly protocol

# optional flags need default values to avoid errors
x=""
a=""
h="true"
L="False"

# define flags
while getopts ":s:T:l:d:x:a:L:" opt; do
  case $opt in
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
    L)
      L=$OPTARG # is there a ligand attached? 
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

# make run folder - assigns UNIX timestamp to name
D=`date +"%s"`
mkdir run${D}
cp "${domains[@]}" ${l} ${x} run${D}/
cd run${D}
mkdir input_scaffold

# prepare scaffold does the following (EC to IC):
# 1) For each input domain, check if it's a dimer (to see where crossing points need to be for design) - hence you need to notify where the ligand is with L
# For now this this means we're limited only to homodimers, if you want to extend it you'll need a clever way of auto-identifying whether it is a dimer (or do everything by hand)
# 2) Removes the linkers based on above arguments (but include them in fasta)
# 3) Prepares topology of input PDBs and writes to all files needed for assembly other than the fragment files

${SCRIPT_DIR}/prepare_scaffold.py -s "${domains[@]}" -d "${d[@]}" -l ${l} ${a} ${x} -L ${L}

mv cst input_scaffold/ 
find -name "*_cut.pdb" -type f -exec mv -t input_scaffold/ {} +
find -name "*fasta" -type f -exec mv -t input_scaffold/ {} +

# rename domains to respect the addition of the _cut suffix
domains_cut=()
for i in "${domains[@]}"; do
    name=$(echo $i | cut -d'.' -f1)
    domains_cut+=("input_scaffold/${name}_cut.pdb")
done

set +o noglob
# finish this, and add script that changes all chain names to A
for s in input_scaffold/*_cut.pdb; do
    python ${SCRIPT_DIR}/renumber_pdb.py ${s}
    tac ${s} | awk '/TER/ {if (f) next; f=1}1' | tac > tmp # remove duplicate TER lines
    mv tmp ${s}
done

printf -v domains_str ' %s' "${domains_cut[@]}"
echo $domains_str

echo "      "
echo ">> Prep completed."
echo ">> Now use the fasta file to create your own fasta fragments (e.g. using https://robetta.bakerlab.org/fragmentsubmit.jsp) before running the mp_assembly on a HPC"

