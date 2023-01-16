#!/bin/bash

# Assess the dimerisation propensity (i.e. energy of binding)
# Q., of the whole system? Or just the domains... likely whole system no?
# use interface analyser

export LD_LIBRARY_PATH="/home/chatzi/Coupl_Ros/ROSETTA.v.3.10_Raj/main/source/build/src/release/linux/5.8/64/x86/gcc/9/default:$LD_LIBRARY_PATH"
R="/home/chatzi/Coupl_Ros/ROSETTA.v.3.10_Raj"
interface="${R}/main/source/bin/InterfaceAnalyzer.linuxgccrelease"
scorer="${R}/main/source/bin/score_jd2.linuxgccrelease"
DATABASE="${R}/main/database/"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# need span file - should be created by relax routine
for i in *.pdb; do
    # save any existing score files
    mv score.sc score_old.sc

    name=$(echo $i | cut -d'.' -f1-3)

    # split chains so interface energy works
    python ${SCRIPT_DIR}/split_chains.py ${i}

    # get both complex energy and interface score 
    $scorer -database ${DATABASE} \
	    -in:membrane true \
	    -score:weights ref2015_memb \
	    -in:file:s ${i} \
	    -mp::setup::spanfiles TM.span \
            -ignore_unrecognized_res false \
	    -ignore_zero_occupancy false \
	    -out:file:score_only ${name}_complex.sc
   
    # franklin2019 doesn't work because of error bad line in file /home/chatzi/Coupl_Ros/ROSETTA.v.3.10_Raj/main/database/scoring/weights/franklin2019.wts:fa_water_to_bilayer 0.5
    # Interface analyser can't use membrane movers...
    $interface -database ${DATABASE} \
	       -score:weights ref2015 \
	       -in:file:s ${i} \
               -ignore_unrecognized_res false \
               -ignore_zero_occupancy false \
               -compute_packstat true \
	       -packstat::oversample 100 \
	       -interface A_B \
	       -out:file:score_only ${name}_interface.sc

               #-score:weights ref2015_memb \
               #-mp::setup::spanfiles TM.span \
               #-in:membrane true \
               #-membrane:Membed_init true \
               #-membrane::Mhbond_depth true \


done
