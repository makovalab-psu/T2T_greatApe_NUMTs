#!/bin/bash

# Arang asks for these NUMT tracks in BigBed format to add to Genome Browser

ASSEMBLY=$1

# Parameters for BigBed tracks.
CHROM_SIZES="chrom_sizes/${ASSEMBLY}.bothHaps.fai"
BB_IN_PRI="numt_set/numts.${ASSEMBLY}.pri.bed"
BB_IN_ALT="numt_set/numts.${ASSEMBLY}.alt.bed"
BB_OUT="../ucsc_*/${ASSEMBLY}.numts.bb"

if [[ $ASSEMBLY == "CHM13" ]]; then CHROM_SIZES="chrom_sizes/${ASSEMBLY}.fai"; BB_IN_ALT=""; fi

# Remove annotations from column 4 and 5.
cat $BB_IN_PRI $BB_IN_ALT | bedtools sort | awk 'BEGIN {OFS="\t"} {$4 = "."; $5 = "0"; print}' | \
   awk 'BEGIN {OFS="\t"} {$4 = "."; $5 = "0"; $6 = "."; print $1, $2, $3, $4, $5, $6}' > temp.bed
# Export as bigBed for the UCSC annotation tracks.
bedToBigBed temp.bed $CHROM_SIZES $BB_OUT


