#!/bin/bash

# change source directories as appropiate
#ROS="/data/rosetta_bin_linux_2020.08.61146_bundle/main/source/build/src/release/linux/3.10/64/x86/gcc/4.8/static"
ROS="/data/rosetta20_glis/source/build/src/debug/linux/5.8/64/x86/gcc/9.3/default"
DATABASE="/data/rosetta_bin_linux_2020.08.61146_bundle/main/database"
R="/data/rosetta20_glis"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
python ${SCRIPT_DIR}/../lib/split_chains.py ${1}

# gen span file
$R/source/bin/mp_span_from_pdb.linuxgccdebug -in:file:s $1
TM=$(echo $1 | rev | cut -d. -f2- | rev)
mv ${TM}.span TM.span

# use relax and idealize alias
idealize="${ROS}/idealize_jd2.default.linuxgccrelease -database ${DATABASE}  -ignore_unrecognized_res false -ignore_zero_occupancy false -overwrite -membrane:Membed_init -membrane::Mhbond_depth -score:weights membrane_highres -spanfile TM.span -in:file:s"
relax="${ROS}/relax.default.linuxgccrelease -database ${DATABASE} -constrain_relax_to_start_coords -relax:ramp_constraints true -relax:bb_move true -nstruct 1 -relax:default_repeats 300 -ignore_unrecognized_res false -ignore_zero_occupancy false -overwrite -membrane:Membed_init -membrane::Mhbond_depth -score:weights membrane_highres -spanfile TM.span -in:file:s"

name=$(echo $1 | cut -d'.' -f1-3)

# first idealize, then relax for a few rounds
$idealize ${name}.pdb
mv ${name}_0001.pdb ${name}_0.pdb

cnt=1
for i in {0..10}; do
    $relax ${name}_${i}.pdb
    mv ${name}_${i}_0001.pdb ${name}_${cnt}.pdb
    cnt=$(($cnt + 1))
done

