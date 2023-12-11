#!/bin/bash
#SBATCH --time 2-12:00:00 -n 1
#SBATCH --array=1-300 # replace with 1-i*j
#SBATCH --cpus-per-task=1
#SBATCH --output slurm.out
#SBATCH --error slurm.err
#SBATCH --partition=serial
#SBATCH --job-name=assemble2

# i = number of repeats you want
# j = number of cluster centres

TM=       # Where are the TMs in your input (likely 2)
R=        # location of Rosetta, e.g. /path/to/Rosetta/main
domain0=  # what is the name of your input file?

output_num=$(echo $domain0 | cut -d'/' -f1 | cut -d'_' -f3)
touch output_scaffold_${output_num}
rm -r output_scaffold_${output_num}
mkdir output_scaffold_${output_num}
O=output_scaffold_${output_num}
I=input_scaffold_${output_num}

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

declare -a commands

for i in {0..9} # how many loops do you want to do?
do
    for j in {0..29} # replace 29 with j (minus 1 to account for counting from 0)
    do
	ij=$(( $i * $n + $j ))
        cmd=$(( $ij + 1 ))

        # now run domain assembly
        # need to resolve static etc. later
        commands[$cmd]="$R/source/bin/mp_domain_assembly.linuxgccrelease \
                        -database $R/database/  \
                        -in:file:fasta ${I}/all.fasta \
                        -in:file:frag3 ${I}/frags.200.3mers \
                        -in:file:frag9 ${I}/frags.200.9mers \
                        -mp:assembly:TM_pose_number ${TM} \
                        -mp:assembly:poses ${domain0} ${I}/c.${j}.0.pdb \
                        -rebuild_disulf false \
                        -detect_disulf false \
                        -constraints::cst_fa_file ${I}/cst \
                        -constraints:cst_fa_weight 1 \
                        -constant_seed \
                        -seed_offset $i \
                        -nstruct 100 \
                        -out:path:all ${O}/ \
                        -out:suffix _${i}_${j} \
                        -out:file:silent out_${i}_${j}.silent \
                        -ignore_zero_occupancy false \
                        > ${O}/log_${i}_${j} \
                        2> ${O}/err_${i}_${j}"
    done
done

bash -c "${commands[${SLURM_ARRAY_TASK_ID}]}"

