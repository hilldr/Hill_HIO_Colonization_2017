#!/bin/bash
## Count number of unique lines in each cell count index file
FILES=dapi_ki67_edu_counts/*
echo "count  file" > cell_counts.tsv
for f in $FILES
do
    wc -l $f >> cell_counts.tsv
done	 
#sed -i '$ d' cell_counts.tsv
tail cell_counts.tsv
