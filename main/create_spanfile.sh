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
