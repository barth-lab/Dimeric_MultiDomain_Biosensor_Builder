#! /usr/bin/bash
R=/data/rosetta20_glis
# since we can't pick up weights directly in the score file, be very penalising in energy
$R/source/bin/score_jd2.linuxgccdebug -in:file:silent output_scaffold/combined.silent -membrane:Membed_init -membrane::Mhbond_depth -score:weights membrane_highres -spanfile input_scaffold/TM.span -no_nstruct_label -constraints true -constraints::cst_fa_file input_scaffold/cst -constraints:cst_fa_weight 100 -viol true  
