#!/bin/bash
set -e
# check for violations of constraints in the output silent files
# extract only those solutions that don't violate any constraints

cat out*silent > combined.silent

# We should always take non-violating structures
grep SCORE combined.silent > tmp.tag
# less than 1 for constraint for a tiny amount of leeway/machine error
awk '$17<1' tmp.tag | awk '{print $NF}' > score.tag
rm tmp.tag

# now extract from silent file, the lines we want given by tag - this includes mem energies we need
# first copy first three lines
head -n 3 combined.silent >> combined_clean.silent
grep -F -f score.tag combined.silent >> combined_clean.silent

# now filter to extract the top 10 % of structures
grep SCORE combined_clean.silent > clean_scores.sc
awk '!/ref/' clean_scores.sc > temp && mv temp clean_scores.sc

# get length of file (to know how many approx 10% represent)
length=$(wc -l clean_scores.sc | awk '{print $1}')
slice_float=$(bc <<< "${length} * 0.1")
slice=${slice_float%.*}

# now get lowest 10% of structures (from column 2 in total score) in new file
sort -nk 2 clean_scores.sc | head -${slice} > filtered_scores.sc
awk '{print $NF}' filtered_scores.sc > filtered_scores.tag

head -n 3 combined_clean.silent >> filtered.silent
grep -F -f filtered_scores.tag combined_clean.silent >> filtered.silent

# cleanup
rm combined_clean.silent combined.silent clean_scores.sc score.tag filtered_scores.tag
