#!/bin/bash
# prep for extraction of build loops - take top 10% of structures again and recluster

cat out*silent > combined.silent

grep SCORE combined.silent > scores.sc

# now filter to extract the top 10 % of structures
awk '!/score/' scores.sc > temp && mv temp top_scores.sc

# get length of file (to know how many approx 10% represent)
length=$(wc -l top_scores.sc | awk '{print $1}')
slice_float=$(bc <<< "${length} * 0.2")
slice=${slice_float%.*}

# now get lowest 10% of structures (from column 2 in total score) in new file
sort -nk 2 top_scores.sc | head -${slice} > filtered_scores.sc
awk '{print $NF}' filtered_scores.sc > filtered_scores.tag

head -n 3 combined.silent >> filtered.silent
while read p; do
    grep $p combined.silent >> filtered.silent
done < filtered_scores.tag

bash /data/domain_construction/domain_assembly_constraints/main/determine_cluster_radius.sh -R /data/rosetta20_glis/ -S combined.silent #-S filtered.silent
