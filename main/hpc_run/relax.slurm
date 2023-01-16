#!/bin/bash
#SBATCH --time 71:59:59 -n 1
#SBATCH --output slurm.out
#SBATCH --error slurm.err
#SBATCH --job-name=Relax
#SBATCH --partition=serial
#SBATCH --cpus-per-task=1
#SBATCH --array=1-12

R=/home/dclw/rosetta20_glis/
n=12 # double the number of structures to loop through (+1)

I=input_relax
O=output_relax
mkdir ${O}

# input structures need to be appropiately names:
# c.0.0, c.1.0, c.2.0, c.3.0 etc.

# same as n above
for i in {1..12}
do
    printf -v sj "%1d" $(((i-1)))

    commands[$i]="$R/source/bin/relax.linuxgccdebug \
    -score:weights membrane_highres \
    -use_input_sc -ex1 -ex2 -extrachi_cutoff 1 \
    -membrane:Membed_init -membrane:Mhbond_depth \
    -database $R/database/  \
    -in:file:s $I/c.${sj}.0.pdb \
    -native $I/c.${sj}.0.pdb \
    -relax::range::angle_max 0.5 \
    -mp:setup:spanfiles $I/TM.span \
    -ignore_unrecognized_res false -ignore_zero_occupancy false \
    -in:membrane \
    -relax:default_repeats 15 \
    -constant_seed \
    -seed_offset $i \
    -nstruct 20 \
    -out:path:all $O \
    -out:suffix _${sj} \
    -out:file:silent out_${sj}.silent \
    > $O/log_${sj} \
    2> $O/err_${sj}"
done

bash -c "${commands[${SLURM_ARRAY_TASK_ID}]}"
    

