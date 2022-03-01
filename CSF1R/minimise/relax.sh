#!/bin/bash

ROS="/data/rosetta_bin_linux_2020.08.61146_bundle/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static"
DATABASE="/data/rosetta_bin_linux_2020.08.61146_bundle/main/database"

# use relax and idealize alias
#idealize="${ROS}/idealize_jd2.static.linuxgccrelease -database ${DATABASE} -output_torsions -suffix _ideal -no_nstruct_label -s"
# Prerelax structure
#relax="${ROS}/relax.static.linuxgccrelease -database ${DATABASE} -constrain_relax_to_start_coords -relax:ramp_constraints false -relax:coord_constrain_sidechains -suffix _relaxed -no_nstruct_label -s"

idealize="${ROS}/idealize_jd2.static.linuxgccrelease -database ${DATABASE}  -ignore_unrecognized_res false -ignore_zero_occupancy false -overwrite -in:file:s"

relax="${ROS}/relax.static.linuxgccrelease -database ${DATABASE} -constrain_relax_to_start_coords -relax:ramp_constraints true -relax:bb_move true -nstruct 1 -relax:default_repeats 300 -ignore_unrecognized_res false -ignore_zero_occupancy false -overwrite -in:file:s"

name=$(echo $1 | cut -d'.' -f1 | cut -d'_' -f 1-2)

cnt=1

for i in {0..10}; do
    $idealize ${name}_${i}.pdb
    mv ${name}_${i}_0001.pdb ${name}_${i}_ideal.pdb
    $relax ${name}_${i}_ideal.pdb
    mv ${name}_${i}_ideal_0001.pdb ${name}_${cnt}.pdb
    rm *_ideal*.pdb
    cnt=$(($cnt + 1))
done

