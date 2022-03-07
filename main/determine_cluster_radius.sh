#!/bin/bash

while getopts ":R:" opt; do
  case $opt in
    R)
      R=$OPTARG # rosetta location
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

round() {
  printf "%.${2}f" "${1}"
}

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# first create a giant silent file from the output of the assembly stage
#cat output_scaffold/out*.silent > output_scaffold/combined.silent # do we need to remove non unique solutions - will effect clustering

# run cluster on the silent output for 10 seconds before crashing and grepping radius from output
$SCRIPT_DIR/cluster.sh $R output_scaffold/combined_clean.silent -1 output_scaffold   &
PID=$(echo $!)
sleep 180
kill ${PID}
big_radius=$(grep --text radius output_scaffold/clout)
echo "big radius size is $big_radius"

# now run in parallel different cluster radius based on big radius - "50% to 80% of the radius used by "-1"."
size_var=(0.5 0.533 0.567 0.6 0.633 0.667 0.7 0.733 0.767 0.8)

for s in ${size_var[@]}; do
    rad=$(bc <<< "${big_radius} * ${s}")
    num=$(round $rad 2) # round to 2 dp
    mkdir output_scaffold/${num}
done


# Then, I run in parallel about 10 clusterings with different radii, from 50% to 80% of the radius used by "-1". 
#I then count how many clusters I get and how many structures there are in the biggest cluster for each radius, 
# and choose the one which seems the most reasonable, i.e. at least 20 clusters, between 10% and 20% of all structures in the largest cluster, 
# at least 1-5% in the smallest cluster, ... Basically, it's educated guessing.

# now run 10 parallel sessions of clustering to see what to ultimatly choose
# at least 20 clusters... so arbitrary - this is something that really needs resolving (e.g. with DBSCAN)
# this is something to implement if you get the funding
