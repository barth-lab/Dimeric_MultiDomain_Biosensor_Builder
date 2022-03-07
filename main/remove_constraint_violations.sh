#!/bin/bash
# check for violations of constraints in the output silent files
# extract only those solutions that don't violate any constraints

cat output_scaffold/out*silent > output_scaffold/combined.silent

grep SCORE output_scaffold/combined.silent > output_scaffold/tmp.tag
# less than 1 for constraint for a small amount of leeway
awk '$17<1' output_scaffold/tmp.tag | awk '{print $NF}' >> output_scaffold/score.tag
rm output_scaffold/tmp.tag

# now extract from silent file, the lines we want given by tag - this includes mem energies we need
# first copy first three lines
head -n 3 output_scaffold/combined.silent >> output_scaffold/combined_clean.silent
while read p; do
    grep $p output_scaffold/combined.silent >> output_scaffold/combined_clean.silent
done < output_scaffold/score.tag
