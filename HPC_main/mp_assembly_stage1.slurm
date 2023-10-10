#!/bin/bash
#SBATCH --time 2-12:00:00 -n 1
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --output slurm.out
#SBATCH --error slurm.err
#SBATCH --job-name=MpAssem1
#SBATCH --partition=serial

O=output_scaffold
mkdir $O

TM=
R=/location/to/rosetta
domains=

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

declare -a commands

for i in {1..100}
do
    printf -v s "%03d" $i

    # now run domain assembly
    # need to resolve static etc. later
    commands[$i]="$R/source/bin/mp_domain_assembly.linuxgccrelease \
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
	-seed_offset $i \
        -nstruct 100 \
        -out:path:all output_scaffold/ \
	-out:suffix _$s \
        -out:file:silent out_$s.silent \
        -ignore_zero_occupancy false \
        > output_scaffold/log_$s \
        2> output_scaffold/err_$s"
done

bash -c "${commands[${SLURM_ARRAY_TASK_ID}]}"
