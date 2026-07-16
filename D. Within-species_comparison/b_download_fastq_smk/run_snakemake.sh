#!/bin/bash

set -eu


# Create log directory if it doesn't exist
mkdir -p logs/slurm

# Initialize micromamba
export MAMBA_ROOT_PREFIX="/storage/home/ejt89"
eval "$($MAMBA_ROOT_PREFIX/bin/micromamba shell hook -s bash)"

# Activate the specific environment
micromamba activate download_fastq




# Run 20 jobs at a time on SLURM.
snakemake -j 20 \
--executor slurm \
--use-conda \
--default-resources mem_mb=16000 runtime=240 \
--output-cache logs/slurm

