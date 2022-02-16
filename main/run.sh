#!/bin/bash

bash /data/domain_construction/domain_assembly_constraints/main/protocol.sh -R /data/rosetta/ -T 4 -s "D1.pdb D2.pdb D3.pdb D4.pdb" -l linkers.txt -d "3 4" -a 3 -x extra_linkers_test.txt
