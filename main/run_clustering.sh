#!/bin/bash

# optional flags need default values to avoid errors
C="-1"
p=0
while getopts ":R:S:C:p:" opt; do
  case $opt in
    R)
      R=$OPTARG # rosetta location
      ;;
    S)
      S=$OPTARG # silent file input
      ;;
    C)
      C=$OPTARG # cluster radius size
      ;;
    p)
      p=$OPTARG
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

$R/source/bin/cluster.linuxgccrelease \
        -database $R/database/  \
        -silent_read_through_errors \
        -in:file:fullatom \
        -in:file:silent $S \
        -cluster:radius $C \
        -run:shuffle \
        -cluster:population_weight 1 \
        -cluster:sort_groups_by_energy \
        -cluster:write_centers \
        > clout \
        2> clerr

if [ "$p" -ne 0 ]
then
    big_radius=$(grep --text radius clout)
    echo "     "
    echo ">> Radius size selected by Rosetta is $big_radius"

    size_var=(0.5 0.533 0.567 0.6 0.633 0.667 0.7 0.733 0.767 0.8 0.833 0.867)
    final_radii=""
    for s in ${size_var[@]}; do
        rad=$(bc <<< "${big_radius} * ${s}")
	num=$(round $rad 2) # round to 2 dp
	final_radii=$(echo "$final_radii $num")
    done

    echo " "
    echo ">> Please rerun the clustering with the following radii sizes: $final_radii" 
fi
