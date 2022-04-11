#!/bin/bash

export LD_LIBRARY_PATH="/home/chatzi/Coupl_Ros/ROSETTA.v.3.10_Raj/main/source/build/src/release/linux/5.8/64/x86/gcc/9/default:$LD_LIBRARY_PATH"
R="/home/chatzi/Coupl_Ros/ROSETTA.v.3.10_Raj"
DCCM="${R}/main/source/bin/dccm_wensel.mpi.linuxgccrelease"

for i in *.pdb; do
    name=$(echo $i | cut -d'.' -f1-3)
    $DCCM -s ${i} > dccm_tmp.txt  
    sed -n '/DCCM_start/,/DCCM_end/{/DCCM_start/b;/DCCM_end/b;p}' dccm_tmp.txt > dccm_matrix.txt
    mv dccm_matrix.txt ${name}_matrix.dccm
    rm dccm_tmp.txt
done
