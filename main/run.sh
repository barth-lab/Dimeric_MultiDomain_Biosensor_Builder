#!/bin/bash

# CD28 test
bash /data/domain_construction/domain_assembly_constraints/main/protocol.sh -R /data/rosetta20_andreas -T 4 -s "D1.pdb D2.pdb D3.pdb D4.pdb" -l linkers.txt -d "3 4" -a 3 -x extra_linkers_test.txt

# VEGFR test
#bash /data/domain_construction/domain_assembly_constraints/main/protocol.sh -R /data/rosetta20_andreas -T 5 -s "D23_dim.pdb  D4.pdb  D5_dim.pdb  D6.pdb  D7TM_dim.pdb" -l linkers.txt -d "1 3 5" -x extra_linkers.txt