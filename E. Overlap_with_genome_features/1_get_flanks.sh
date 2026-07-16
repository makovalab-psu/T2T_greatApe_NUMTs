#!/bin/bash

set -eu
source config.sh


######################################
## Remove hap and hsa from chr name ##
######################################

mkdir -p $NUMT_DIR
for ASS in "${ASSEMBLIES[@]}"; do
	OG_NUMTS=${OG_NUMT_DIR}/numts.${ASS}.pri.bed
	cat $OG_NUMTS | awk '{OFS="\t"}{gsub(/_.*/, "", $1); print $1,$2,$3,$4,$7,$6}' > ${NUMT_DIR}/${ASS}.bed
done


#################
## NUMT flanks ##
#################

# FLANK_SIZE=50

for ASS in "${ASSEMBLIES[@]}"; do
	NUMTS=${NUMT_DIR}/${ASS}.bed
	bedtools flank -i ${NUMTS} -g ${GENOME_DIR}/${ASS}*fai -b $FLANK_SIZE | bedtools sort | bedtools merge -c 4,5,6 -o distinct > ${OUTPUT_DIR}/flanks.${ASS}.bed ; 
done


#################
## NUMT flanks ##
#################

for ASS in "${ASSEMBLIES[@]}"; do
	NUMTS=${NUMT_DIR}/${ASS}.bed
  	TOTAL=$( bedtools merge -i $NUMTS | awk '{total += $3 - $2} END {print total}' ); 
 	echo "${ASS} ${TOTAL}" | sed 's/ /\t/g' ; 
done > ${NUMT_DIR}/totalbp_numts.txt

for ASS in "${ASSEMBLIES[@]}"; do
	NUMTS=${NUMT_DIR}/${ASS}.bed
  	TOTAL=$( bedtools merge -i $NUMTS | wc -l | awk '{OFS="\t"}{$1=$1}{print}' ); 
  	echo "${ASS} ${TOTAL}" | sed 's/ /\t/g' ; 
done > ${NUMT_DIR}/counts_numts.txt

