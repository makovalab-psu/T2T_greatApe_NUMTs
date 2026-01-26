#!/bin/bash


# Parameters for overlapping with genomic features
FLANK_SIZE=50           # NUMT flanks size in base pairs
INTERSECTION_RATIO=0.5  #0.00001   # At least 80% of the NUMT must overlap with annotation


# Directories for NUMTs and genomes
OG_NUMT_DIR="../../numts_from_laptop_blast/numt_set"
NUMT_DIR="numts"
mkdir -p numts
GENOME_DIR="pri_species_chr"

# List of assemblies to process
ASSEMBLIES=(
  "mGorGor1"
  "mPanPan1"
  "mPanTro3"
  "mPonAbe1"
  "mPonPyg2"
  # "mSymSyn1"
  "CHM13"
)

# Directories for genome features
FEATURE_DIR="../../../bin/F_numts_and_genes/functionalAnnotations"
# FEATURES=( "intron" "exon" "gene_protCode" "gene_nonprotCode" "utr3" "utr5" "cds" "rna" "enhancer" "promoter" "chrG_cpgislands" "ngnr" "repeats" "te" "nonTe" "nf_te" "simple" "satellite" )
FEATURES=( "ngnr" "repeats" "te" "nonTe" "nf_te" "simple" "satellite" )

# BLAT output directories
OUTPUT_DIR="./output_${FLANK_SIZE}bp"
OVERLAPS_DIR="${OUTPUT_DIR}/ir_${INTERSECTION_RATIO}"
LOG_DIR="./log"

# Ensure required directories exist
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$OVERLAPS_DIR"

# Logging configuration
echo "Processing started. Logs will be stored in $LOG_DIR."
echo "Output files will be stored in $OUTPUT_DIR."
