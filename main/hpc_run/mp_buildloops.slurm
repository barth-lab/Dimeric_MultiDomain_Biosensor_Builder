#!/bin/bash
#SBATCH --time 2-12:00:00 -n 1
#SBATCH --cpus-per-task=1
#SBATCH --output slurm.out
#SBATCH --error slurm.err
#SBATCH --job-name=BuildLoops
#SBATCH --partition=serial
#SBATCH --array=1-361
# batch scripts don't allow this to be precalculated sadly

R=/home/dclw/rosetta20_glis/
n=17 # number of cluster center PDBs to use (technically +1, but start at c.0.0.pdb)

# number of array jobs is n +1 * 20

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

O=output_loop
I=input_loop

declare -a commands

for i in {0..19}
do
    for j in {0..17}; do

        ij=$(( $i * $n + $j ))
        cmd=$(( $ij + 1 ))

        commands[$cmd]="$R/source/bin/loopmodel.linuxgccdebug \
                -database $R/database/  \
                -in:file:s $I/c.${j}.0.pdb \
                -loop_file $I/loopfile \
                -loops:frag_sizes 9 3 \
                -loops:frag_files input_scaffold/frags.200.9mers input_scaffold/frags.200.3mers \
                -loops:remodel quick_ccd \
                -loops:refine refine_ccd \
                -loops:relax no \
                -nstruct 10 \
                -run:use_time_as_seed \
                -run:seed_offset $cmd \
                -out:suffix _${j}_${i} \
	        -out:file:silent out_${i}_${j}.silent \
                -out:path:all $O \
                -ignore_zero_occupancy false \
                > $O/log \
                2> $O/err"
    done
done

bash -c "${commands[${SLURM_ARRAY_TASK_ID}]}"
