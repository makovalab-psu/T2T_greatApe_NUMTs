#!/bin/bash

set -eu

# '''
# Reproduce the NUMT discovery pipeline using BLASTn.
# '''

# Arguments:
#ASS=$1
#HAPLOTYPE=$2 # "pri" or "alt"


for ASS in "mGorGor1" "mPanPan1" "mPanTro3" "mPonAbe1" "mPonPyg2" "mSymSyn1" "CHM13" ; do
 echo "-$ASS"
 for HAPLOTYPE in "pri" "alt" ; do
  echo "--$HAPLOTYPE"

# Parameters:
REFDIR="../../data/refs"
NUCL="${REFDIR}/${ASS}.${HAPLOTYPE}.*.fasta"
if [[ $HAPLOTYPE == 'alt' ]]; then NUCL="${REFDIR}/alternative/${ASS}.${HAPLOTYPE}.cur.*.fasta"; fi
if [[ $ASS == 'CHM13' ]]; then NUCL="${REFDIR}/*${ASS}*.fna"; fi
if [[ $HAPLOTYPE == 'alt' && $ASS == 'CHM13' ]]; then echo "# There is no alternate haplotype for haploid CHM13" ; continue; fi  
MT="${REFDIR}/mt_genomes/doubled_genome/${ASS}*.doubled.fasta"
OUTDIR="blast_results"
mkdir -p blast_results
NUMT_TAB=${OUTDIR}/blast.${ASS}.${HAPLOTYPE}.tab
NUMT_BED=${OUTDIR}/blast.${ASS}.${HAPLOTYPE}.bed

# # If the DB does not have an index file.
# cd $REFDIR
# makeblastdb -in $NUCL -dbtype nucl -out $NUCL

# # Double the MT genome (to account for ref breakpoint).
# cd $OUTDIR
# python double_fasta.py $MT
# MT=${OUTDIR}/${ASS}*.doubled.fasta

# Align mtDNA to T2T nuclear assembly.
blastn  -db $NUCL -query $MT \
        -outfmt '7' -task blastn -evalue 0.0001 \
        -gapopen 5 -gapextend 2 -penalty -3 -reward 2 \
        > $NUMT_TAB

# Process blastn output file into BED format.
# # The BED coordinates were sorted then merged overlaps (accounts for doubled genome).
python process_blast.py $NUMT_TAB | bedtools sort | bedtools merge -c 4,5,6 -o distinct > $NUMT_BED #candidate NUMTs

# Get coordinates for the 500 bp flanks alongside NUMTs.
bedtools flank -b 500 -i $NUMT_BED -g ${NUCL}.fai > ${OUTDIR}/flanks.${ASS}.${HAPLOTYPE}.bed

# Retrieve flank sequences.
bedtools getfasta -name -fi $NUCL -bed ${OUTDIR}/flanks.${ASS}.${HAPLOTYPE}.bed -fo ${OUTDIR}/flanks.${ASS}.${HAPLOTYPE}.fasta

# Concatenate flank sequences.
python concat_fasta.py ${OUTDIR}/flanks.${ASS}.${HAPLOTYPE}.fasta

# BLAT Sorang NUMTs to Sorang genome.
Q=$ASS
DB=$ASS
QUERY=${OUTDIR}/concat.flanks.${Q}.${HAPLOTYPE}.fasta
DATABASE=${REFDIR}/${DB}.${HAPLOTYPE}.cur.*.fasta
#time blat -minScore=650 $DATABASE $QUERY ${OUTDIR}/matches.${Q}.to.${DB}.${HAPLOTYPE}.psl

 done

done


# Copy NUMTs to a common directory.
mkdir -p numt_set
# Add a unique identifier to each NUMT.
for ASS in "mGorGor1" "mPanPan1" "mPanTro3" "mPonAbe1" "mPonPyg2" "mSymSyn1" ; do for HAPLOTYPE in "alt" "pri"; do awk -v prefix="${ASS}_${HAPLOTYPE}_N" 'BEGIN { OFS = "\t" } { print $0, prefix i++ }' blast_results/blast.$ASS.${HAPLOTYPE}.bed > numt_set/numts.$ASS.${HAPLOTYPE}.bed  ; done; done; 
awk -v prefix="CHM13_pri_N" 'BEGIN { OFS = "\t" } { print $0, prefix i++ }' blast_results/blast.CHM13.pri.bed > numt_set/numts.CHM13.pri.bed


