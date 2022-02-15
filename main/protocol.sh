#!/bin/bash

# author: Lucas Rudden

# Here I want to write out the entire protocol of the modified approach to domain assembly
# Assembly procotol is essentially domain assemble (with constraints), clustered and then taken further for assembly rounds
# Construction always needs to occur down (i.e. towards cytoplasm)

O=$(pwd)

# make run folder
D=`date +"%s"`
mkdir run${D}
cd run${D}

# Input domains - make this generalisable later with number of input arguments
D1=${O}/$1
D2=${O}/$2
D3=${O}/$3
D4=${O}/$4

# for each input domain, check if it's a dimer (to see where crossing points need to be for design) - note, this need to avoid checking for ligand!
# FOr now this can be done by seeing if its a homodimer (how to do heterodimers later?)
# Then we need to remove the linkers (but include within fasta)
# then we need to build the fragment database
# THen we run domain assembly

# OTHER NEEDED INPUT FILES:
# 1) fragment files (for domain assembly and linker modelling) - can be built with fasta file created by prepare_scaffold
# 2) fasta file of domains and, eventually, chimera
# 3) Span file for the TM domain

# two protocols - one to run on cluster, one without?

# Starting scaffold needs the start and end points included (i.e. EC head to CT region) to build inial scaffold. This is based on the input domains
python prepare_scaffold.py -s D1.pdb D2.pdb D3.pdb D4.pdb -d 3 4 -a 3 -l linkers.txt -x extra_linkers_test.txt 

# need to cut the linker regions during assembly as they'll be readded later - however we don't want to apply this to D3 as this is unstructured already

# run domain assembly with constraints - building from top to bottom (EC to CT)

# cluster

# run linker design

# cluster?

# select from miminum energy, and repear
