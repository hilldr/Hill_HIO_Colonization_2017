#!/bin/bash
## AUTHOR: David R. Hill
## Count number of unique lines in each cell count index file
FILES=../data/figure4_/dapi_ki67_edu_counts/*
echo "count  file" > ../data/figure4_/cell_counts.tsv
for f in $FILES
do
    wc -l $f >> ../data/figure4_/cell_counts.tsv
done	 
