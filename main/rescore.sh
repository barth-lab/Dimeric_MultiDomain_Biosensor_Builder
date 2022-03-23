#! /usr/bin/bash
R=/data/rosetta20_glis
# make spanfile
$R/source/bin/mp_span_from_pdb.linuxgccdebug -in:file:s $1
TM=$(echo $1 | rev | cut -d. -f2- | rev)
mv ${TM}.span TM.span
# since we can't pick up weights directly in the score file, be very penalising in energy
# Silent
#$R/source/bin/score_jd2.linuxgccdebug -in:file:silent $1 -membrane:Membed_init -membrane::Mhbond_depth -score:weights membrane_highres -spanfile TM.span -no_nstruct_label -constraints true -constraints::cst_fa_file cst -constraints:cst_fa_weight 1 -viol true  
# PDB
$R/source/bin/score_jd2.linuxgccdebug -in:file:s $1 -membrane:Membed_init -membrane::Mhbond_depth -score:weights membrane_highres -spanfile TM.span -no_nstruct_label
