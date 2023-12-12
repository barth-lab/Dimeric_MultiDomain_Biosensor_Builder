#!/bin/bash

S=2
N=3
# define flags
while getopts ":R:S:T:N:D:" opt; do
  case $opt in
    R)
      R=$OPTARG # Location of Rosetta
      ;;
    S)
      S=$OPTARG # stage of assembly
      ;;
    T)
      TM=$OPTARG # which domain contains the TM?
      ;;
    N)
      N=$OPTARG # number of structures to generate (default=3)
      ;;
    D)
      domain0=$OPTARG # domain to add in
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

touch output_scaffold_${S}
rm -r output_scaffold_${S}
mkdir output_scaffold_${S}
O=output_scaffold_${S}
I=input_scaffold_${S}

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

for i in ${I}/c.*.0.pdb; do
    idx=$(echo $i | cut -d'/' -f2 | cut -d'.' -f2)
    $R/source/bin/mp_domain_assembly.linuxgccrelease \
              -database $R/database/  \
              -in:file:fasta ${I}/all.fasta \
              -in:file:frag3 ${I}/frags.200.3mers \
              -in:file:frag9 ${I}/frags.200.9mers \
              -mp:assembly:TM_pose_number ${TM} \
              -mp:assembly:poses ${domain0} ${i} \
              -rebuild_disulf false \
              -detect_disulf false \
              -constraints::cst_fa_file ${I}/cst \
              -constraints:cst_fa_weight 1 \
              -constant_seed \
              -nstruct $N \
              -out:path:all ${O}/ \
              -out:file:silent out${idx}.silent \
              -ignore_zero_occupancy false \
              > ${O}/log${idx} \
              2> ${O}/err${idx}
done
