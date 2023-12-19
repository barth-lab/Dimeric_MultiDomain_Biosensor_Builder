#!/bin/bash

N=6
# define flags
while getopts ":R:T:d:N:" opt; do
  case $opt in
    T)
      TM=$OPTARG # which input domain corresponds to the TM region
      ;;
    d)
      # allow importing of multiple domains https://unix.stackexchange.com/questions/164259/provide-two-arguments-to-one-option-using-getopts 
      set -f # disable glob
      IFS=' ' # split on space characters
      domains=($OPTARG) ;; # domain pdbs
    R)
      R=$OPTARG # Location of Rosetta
      ;;
    N)
      N=$OPTARG # number of output pdbs you're asking for
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

O=output_scaffold
touch $O
rm -r $O
mkdir $O

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

$R/source/bin/mp_domain_assembly.linuxgccrelease \
	      -database $R/database/  \
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
	      -nstruct $N \
	      -out:path:all output_scaffold/ \
	      -ignore_zero_occupancy false \
	      -out:file:silent out.silent \
              > output_scaffold/log \
              2> output_scaffold/err
