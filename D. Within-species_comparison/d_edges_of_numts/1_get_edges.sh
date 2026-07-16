#!/bin/bash

set -eu

REF_DIR="/storage/group/kdm16/default/shared/T2Tv2.assemblies/NCBI_RefSeq"
PONPYG="PonPyg/data/GCF_028885625.2/GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna.gz" 

# Make 100 bp windows at NUMT start and end.
awk 'BEGIN {OFS="\t"} {print $1, $2-50, $2+50, $4, $5, $6, $7, $8}' numts.mPonPyg2.pri.bed > centered_starts.bed
awk 'BEGIN {OFS="\t"} {print $1, $3-50, $3+50, $4, $5, $6, $7, $8}' numts.mPonPyg2.pri.bed > centered_ends.bed

# Get sequence.
bedtools getfasta -fi $REF_DIR/$PONPYG -bed centered_starts.bed | bgzip > centered_starts.fa.gz
bedtools getfasta -fi $REF_DIR/$PONPYG -bed centered_ends.bed | bgzip > centered_ends.fa.gz

# Align windows to sequencing reads.
srun --nodes=1 --cpus-per-task=4 --mem=16G --time=02:00:00 --pty bash -c \
"minimap2 -a -x sr -t 4 centered_starts.fa.gz $REF_DIR/$PONPYG | samtools view -bh - | samtools sort -o centered_starts.bam -"


