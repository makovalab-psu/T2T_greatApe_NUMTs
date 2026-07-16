#!/bin/bash

set -eu

## Use the full, unclipped graph HPRCv1 genome (ref CHM13) to extract NUMT positions.

mkdir -p extracted_regions

#cat numts.CHM13.pri.bed | sed 's/CHM13\.noNC\://g' | awk '{OFS="\t"}{print $8,$2,$3,$7,$5,$6,$4}' > fixed.numts.CHM13.pri.bed
BED_FILE="fixed.numts.CHM13.pri.bed"

# Loop through each region in the BED file
while read -r chrom start end name score strand info; do
  # Fallback if the 4th column (name) is empty
  if [ -z "$name" ]; then name="${chrom}_${start}_${end}"; fi
  
  echo "Extracting structural locus: ${name}..."
  
  # Extract the specific region based on the exact coordinates
  OMP_NUM_THREADS=2 odgi extract \
    -i hprc-v1.1-mc-chm13.full.og \
    -r "${chrom}:${start}-${end}" \
    -o "extracted_regions/${name}.og"
    
done < $BED_FILE

