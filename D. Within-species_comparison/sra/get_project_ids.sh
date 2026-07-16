#!/bin/bash
set -eu

# Create a text file with one SRR ID per line
cat SRR_list.txt | while read srr; do
    echo -n "$srr "
    #esearch -db sra -query "$srr" | efetch -format docsum | grep -o "Bioproject [0-9]*"
    efetch -db sra -id ${srr} -format docsum | grep 'Bioproject'
done


