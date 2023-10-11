#!/bin/bash

big_radius=$1

round() {
  printf "%.${2}f" "${1}"
}

# now run in parallel different cluster radius based on big radius - "50% to 80% of the radius used by "-1"."
size_var=(0.5 0.533 0.567 0.6 0.633 0.667 0.7 0.733 0.767 0.8 0.833 0.867)

final_radii=""
for s in ${size_var[@]}; do
    rad=$(bc <<< "${big_radius} * ${s}")
    num=$(round $rad 2) # round to 2 dp
    final_radii=$(echo "$final_radii $num")
done

echo " "
echo ">> Please rerun the clustering with the following radii sizes: $final_radii"

