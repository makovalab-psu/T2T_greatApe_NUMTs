#!/bin/bash

set -ue
source config.sh

# FLANK_SIZE=$1
# OUTPUT_DIR="./output_${FLANK_SIZE}bp"
# OVERLAPS_DIR="${OUTPUT_DIR}/ir_${INTERSECTION_RATIO}"

#####################################
# Annotate NUMTs with ________
## Keep only full NUMT overlaps with _____

# INTERSECTION_RATIO=0.8 #at least a 80% of the NUMT overlaps with annotation


# for FEATURE in intron exon gene_protCode gene_nonprotCode utr3 utr5 cds rna enhancer promoter chrG_cpgislands ngnr repeats te simple satellite; do

for FEATURE in "${FEATURES[@]}"; do
 echo "-${FEATURE}"

 # Parameters.
 WORKDIR="${OVERLAPS_DIR}/flanks_and_${FEATURE}s"
 mkdir -p $WORKDIR
 DIR_ANNOT=$FEATURE_DIR


 for ASS in "${ASSEMBLIES[@]}"; do
  echo "--${ASS}"

  # Parameters.
  ANNOTATIONS="${DIR_ANNOT}/${ASS}/${FEATURE}*.bed"
  FAI="${GENOME_DIR}/${ASS}*.fai"
  FLANKS=${OUTPUT_DIR}/flanks.${ASS}.bed

  if [[ $FEATURE == 'ngnr' ]]; then ANNOTATIONS="${DIR_ANNOT}/${ASS}/${FEATURE}*.bed.gz"; fi
  if [[ $FEATURE == 'repeats' ]]; then ANNOTATIONS="${DIR_ANNOT}/${ASS}/repeats_protCode*.bed.gz"; fi

  # # Intersect NUMTs and functional annotations
  # ## (-f 1e-9) Is the minimum overlap (1 bp).
  # ### (-f 1.0) Forces complete overlap of NUMT within a repeat annotation as fraction of a.
  # bedtools merge -i $ANNOTATIONS -c 4,5,6 -o distinct | bedtools intersect -F $INTERSECTION_RATIO -a - -b $FLANKS -loj > ${WORKDIR}/${ASS}.${FEATURE}.numts.bed

  # gzip -f ${WORKDIR}/${ASS}.${FEATURE}.numts.bed

  # # Simplify overlaps file to an index, NUMT presence, NUMT_ID, and feature name.
  # gzcat ${WORKDIR}/${ASS}.${FEATURE}.numts.bed.gz | sed 's/\t\-1/\t./g' | cut -f4,11 | uniq | awk '{OFS="\t"} {i+=1} {print i, ($2 == ".") ? 0 : 1, $2, $1}' > ${WORKDIR}/${ASS}.${FEATURE}.numts.simpl

  # gzip -f ${WORKDIR}/${ASS}.${FEATURE}.numts.simpl

  # python run_bootstrap.py ${WORKDIR}/${ASS}.${FEATURE}.numts.simpl.gz
  gzip -f ${WORKDIR}/${ASS}.${FEATURE}.numts.simpl_bootstrap_frequencies.txt

 done

 # ############################
 # # Join all species into one. Merge overlapping intervals. 
 # cd $WORKDIR
 # cat *[0-9].${FEATURE}.numts.bed | grep ':' | bedtools sort | bedtools merge -c 4,5,6 -o distinct > all.${FEATURE}.numts.bed

 # ###################################################
 # # Compute counts and Total BP for each species and annotation.
 
 # # Total BP.
 # for ASS in "${ASSEMBLIES[@]}"; do
 #  TOTAL=$( bedtools sort -i ${ASS}.${FEATURE}.numts.bed | bedtools merge | awk '{total += $3 - $2} END {print total}' ); 
 #   echo "${ASS} ${TOTAL}" | sed 's/ /\t/g' ; 
 # done > totalbp_numts_and_${FEATURE}.txt

 # # Counts.
 # for ASS in "${ASSEMBLIES[@]}"; do
 #  TOTAL=$( bedtools sort -i ${ASS}.${FEATURE}.numts.bed | bedtools merge | wc -l | awk '{OFS="\t"}{$1=$1}{print}' );
 #   echo "${ASS} ${TOTAL}" | sed 's/ /\t/g' ;
 # done > counts_numts_and_${FEATURE}.txt
 # cd ../../..

done



