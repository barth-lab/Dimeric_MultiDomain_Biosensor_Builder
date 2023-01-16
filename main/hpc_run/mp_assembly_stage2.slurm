#!/bin/bash
#SBATCH --time 2-12:00:00 -n 1
#SBATCH --array=1-281
#SBATCH --cpus-per-task=1
#SBATCH --output slurm.out
#SBATCH --error slurm.err
#SBATCH --partition=serial
#SBATCH --job-name=assemble2

O=output_scaffold_2
mkdir $O

TM=2
R=/home/dclw/rosetta20_glis/
domain0="input_scaffold_2/D6_cut.pdb"
n=28 # how many cluster centers do we have (what is the number of the last, so +1)

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

declare -a commands

for i in {0..9}
do
    for j in {0..28}
    do
	ij=$(( $i * $n + $j ))
        cmd=$(( $ij + 1 ))

        # now run domain assembly
        # need to resolve static etc. later
        commands[$cmd]="$R/source/bin/mp_domain_assembly.linuxgccdebug \
            -database $R/database/  \
            -in:file:fasta input_scaffold_2/all.fasta \
            -in:file:frag3 input_scaffold_2/frags.200.3mers \
            -in:file:frag9 input_scaffold_2/frags.200.9mers \
            -mp:assembly:TM_pose_number ${TM} \
            -mp:assembly:poses ${domain0} input_scaffold_2/c.${j}.0.pdb \
            -rebuild_disulf false \
            -detect_disulf false \
            -constraints::cst_fa_file input_scaffold_2/cst \
            -constraints:cst_fa_weight 1 \
            -constant_seed \
            -seed_offset $i \
            -nstruct 100 \
            -out:path:all output_scaffold_2/ \
            -out:suffix _${i}_${j} \
            -out:file:silent out_${i}_${j}.silent \
            -ignore_zero_occupancy false \
            > output_scaffold_2/log_${i}_${j} \
            2> output_scaffold_2/err_${i}_${j}"
    done
done

bash -c "${commands[${SLURM_ARRAY_TASK_ID}]}"

