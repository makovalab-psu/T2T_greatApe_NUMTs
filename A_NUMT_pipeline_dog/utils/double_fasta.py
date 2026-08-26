#!/usr/bin/env python3
import sys
import os
from Bio import SeqIO

# Doubles sequence entries (useful for circular mt genomes to avoid breakpoint issues)

def double_fasta_entry(input_path):
    input_file = os.path.basename(input_path)
    file_suffix = input_file.split('.')[-1]
    if file_suffix in ['fa', 'fasta']:
        file_prefix = '.'.join(input_file.split('.')[:-1])
        output_file = f"{file_prefix}.doubled.{file_suffix}"
    else:
        raise ValueError(f"Need a FASTA file. {input_path} does not end in .fa or .fasta")

    records = list(SeqIO.parse(input_path, "fasta"))
    with open(os.path.join(os.getcwd(), output_file), "w") as outfile:
        for record in records:
            combined_sequence = str(record.seq) * 2
            combined_id = f"{record.id} (doubled)"
            outfile.write(f">{combined_id}\n{combined_sequence}\n")
    print(f"Doubled FASTA written to {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: double_fasta.py <input_fasta>')
        sys.exit(1)
    double_fasta_entry(sys.argv[1])
