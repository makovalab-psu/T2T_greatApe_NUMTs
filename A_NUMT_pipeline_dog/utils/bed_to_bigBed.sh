#!/bin/bash
set -eu
ASSEMBLY=$1
CHROM_SIZES="chrom_sizes/${ASSEMBLY}.bothHaps.fai"
BB_IN_PRI="numt_set/numts.${ASSEMBLY}.pri.bed"
BB_IN_ALT="numt_set/numts.${ASSEMBLY}.alt.bed"
BB_OUT="${ASSEMBLY}.numts.bb"
if [[ $ASSEMBLY == "CHM13" ]]; then CHROM_SIZES="chrom_sizes/${ASSEMBLY}.fai"; BB_IN_ALT=""; fi
cat $BB_IN_PRI $BB_IN_ALT | bedtools sort | awk 'BEGIN {OFS="\t"} {$4 = "."; $5 = "0"; print}' | awk 'BEGIN {OFS="\t"} {$4 = "."; $5 = "0"; $6 = "."; print $1, $2, $3, $4, $5, $6}' > temp.bed
bedToBigBed temp.bed $CHROM_SIZES $BB_OUT
