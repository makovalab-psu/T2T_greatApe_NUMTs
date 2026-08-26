#!/bin/bash

# Pipeline configuration – edit for your environment
# REFDIR should contain assembly FASTA(s) (e.g., canFam3.fasta) and an mt_genomes/ subdir with mt FASTA(s)
REFDIR="/path/to/refs"   # <-- set this before running
ASSEMBLIES_FILE="assemblies.txt"
MT_SIZES_FILE="mt_sizes.tsv"
# Haplotypes to consider; for most dog assemblies a single haplotype 'pri' is sufficient
HAPLOTYPES=("pri")

# Output directory for intermediate BLAST results
OUTDIR="blast_results"
mkdir -p "$OUTDIR"

# Export for scripts
export REFDIR OUTDIR ASSEMBLIES_FILE MT_SIZES_FILE
