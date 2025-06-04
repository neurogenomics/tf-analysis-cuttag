# Only keep read 1s from alignment files for 5'-shift peak calling
# More information: 
# https://www.biorxiv.org/content/10.1101/2024.10.08.617149v1.full

#!/bin/bash

#### Tweakable stuff ###########################################################
INP_FOLDER="data/linear_dedup"
OUT_FOLDER="data/linear_dedup_read1"

#### Run #######################################################################
mkdir -p $OUT_FOLDER

for file in $INP_FOLDER/*.bam; do
    base=$(basename $file)
    echo "Processing $base"
    bamtools filter -in $file -out $OUT_FOLDER/$base -isFirstMate true
done

echo "All read 1s have been filtered and saved to $OUT_FOLDER"