#!/bin/bash
set -e

# An integration test for ensuring the Dimeric MultiDomain Biosensor Builder runs okay
# Uses the minimum prompts possible in the interest of speed
# The only parameter you should need to provide is your local rosetta build location

while getopts ":R:" opt; do
  case $opt in
    R)
      ROS=$OPTARG # Rosetta build location
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# current script location
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

touch TEST
rm -r TEST
mkdir TEST

cp ${SCRIPT_DIR}/example/* ./TEST
cd ./TEST/

echo ">> Crafting scaffold topology files"

${SCRIPT_DIR}/main/prepare_scaffold.sh -T 4 -s "D13.pdb D45.pdb D6.pdb D7.pdb" -l remove_linkers.txt -d "1 2 4" -a "1 2 3 4" -x add_linkers.txt -L True

runfolder=$(find -name "run*")
cd ./$runfolder
cp ../frags* input_scaffold/

orig=$(pwd)

echo ">> Assembling the initial scaffold"

${SCRIPT_DIR}/main/mp_assembly_stage1.sh -T 4 -R ${ROS} -d "input_scaffold/D13_cut.pdb input_scaffold/D45_cut.pdb input_scaffold/D6_cut.pdb input_scaffold/D7_cut.pdb"

echo ">> Removing any constraint violating files"

cd ./out_scaffold/
${SCRIPT_DIR}/main/remove_constraint_violations.sh

echo ">> Clustering"
${SCRIPT_DIR}/main/run_clustering.sh -R ${ROS} -S filtered.silent

mkdir 0.0
mv c.*.pdb 0.0/
cd 0.0
cnt=0
for i in *.pdb; do
    mv $i c.${cnt}.0.pdb
    ((cnt = cnt + 1))
done

cd ${orig}
echo ">> Setting up the scaffold for stage 2 by preparing for the missing empty domain"
${SCRIPT_DIR}/main/assemble_domains.sh -S 2 -d "1 2 4" -s "D6.pdb c.0.0.pdb" -C 0.0 -p 3

echo ">> Constructing the next scaffold"
${SCRIPT_DIR}/main/mp_assemble_stage2.sh -R ${ROS} -S 2 -T 2 -D "input_scaffold_2/D6_cut.pdb"

echo ">> Reclustering"
cd ./output_scaffold_2
${SCRIPT_DIR}/main/remove_constraint_violations.sh
${SCRIPT_DIR}/main/run_clustering.sh -R ${ROS} -S filtered.silent

mkdir 0.0
mv c.*.pdb 0.0/
cd 0.0
cnt=0
for i in *.pdb; do
    mv $i c.${cnt}.0.pdb
    ((cnt = cnt + 1))
done

cd ${orig}

echo ">> Preparing to rebuild the final loops"
${SCRIPT_DIR}/main/prepare_linkers.sh -C 0.0 -d "1 2 4" -o "3 6 7 8 4 5 1 2" -l 3 -S 2

echo ">> Building loops!"
${SCRIPT_DIR}/main/build_linkers.sh -R ${ROS}

echo ">> Integration test complete, check ${orig}/output_loop for final models"

cd ${SCRIPT_DIR}
