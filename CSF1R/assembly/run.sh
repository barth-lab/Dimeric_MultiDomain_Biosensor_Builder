#!/bin/bash

bash /data/domain_construction/domain_assembly_constraints/main/protocol.sh -R /data/rosetta20_glis -T 2 -s "D1.pdb D2.pdb" -l linkers.txt -d "1 2" -a "1 2" -x extra_linkers.txt

