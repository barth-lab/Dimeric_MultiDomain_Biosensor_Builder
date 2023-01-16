#!/bin/bash

#Remove structures with far placed atoms (Rosetta failure state)

if [[ -z $1 ]]
then
    echo "Usage: $0 input.pdb"
    echo "  task: Remove a PDB if it has loop atoms placed very far away"
    exit
fi

farAwayWaters=0
while read -r line
do
    if [[ "${line:0:4}" == 'ATOM' && $(echo "${line:29:9} > 900.0" | bc) == 1 ]]
    then
        farAwayWaters=$(echo "${farAwayWaters} + 1" | bc)
    fi
done < ${1}

if [[ $farAwayWaters > 1 ]]
then
    echo "Found ${farAwayWaters} atoms that are far away in ${1}"
    rm ${1}
fi

