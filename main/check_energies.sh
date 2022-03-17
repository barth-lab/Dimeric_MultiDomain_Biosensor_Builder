#!/bin/bash
# check if scores are less than certain value and pipe to file

# less than 1 for constraint for a small amount of leeway
awk '$2<0' score.sc | awk '{print $NF}' >> low_score.tag

head -n 3 output_scaffold/combined.silent >> output_scaffold/combined_lowE.silent
while read p; do
    grep $p output_scaffold/combined.silent >> output_scaffold/combined_lowE.silent
done < low_score.tag
