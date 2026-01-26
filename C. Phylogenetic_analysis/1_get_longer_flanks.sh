#!/bin/bash
set -ue

#ASS=$1

FLANK_SIZE=500
for ASS in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 CHM13; do 

GENOME_FILE=../../data/refs/${ASS}.pri.cur.*.fasta
if [[ $ASS == "CHM13" ]]; then GENOME_FILE=../../data/*/GCF*CHM13*.fna ; fi

## FLANKS ###
# Keep just the flanks alongside NUMTs.
bedtools flank -b $FLANK_SIZE -i ../numts/merged.blast.rawCoords.${ASS}.bed -g ${GENOME_FILE}.fai > flanks${FLANK_SIZE}.${ASS}.bed
# Retrieve flank sequence.
bedtools getfasta -name -fi $GENOME_FILE -bed flanks${FLANK_SIZE}.${ASS}.bed -fo flanks${FLANK_SIZE}.${ASS}.fasta

done