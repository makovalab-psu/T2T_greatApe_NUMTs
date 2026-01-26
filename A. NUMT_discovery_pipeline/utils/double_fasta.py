
import sys
import os
from Bio import SeqIO

def double_fasta_entry(input_path):
    # Determine the output file name
    input_file = input_path.split('/')[-1]
    file_suffix = input_file.split('.')[-1]
    if file_suffix in ['fa', 'fasta']:
        file_prefix = '.'.join(input_file.split('.')[:-1])
        output_file = f"{file_prefix}.doubled.{file_suffix}"
    else:
        raise ValueError(f"### Need to input a FASTA file. {input_path} does not end in .fa or .fasta ###")
    
    # Parse the input file for multiple entries
    with open(input_path, "r") as infile:
        records = list(SeqIO.parse(infile, "fasta"))
    
    # Write the doubled sequences to output file
    with open(os.path.join(os.getcwd(), output_file), "w") as outfile:
        for record in records:
            combined_sequence = str(record.seq) * 2  # Double the sequence
            combined_id = f"{record.id} (doubled)"   # Update the ID
            outfile.write(f">{combined_id}\n{combined_sequence}\n")
    
    print(f"### Doubled FASTA written to {output_file} ###")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python concat_fasta.py <input_file>")
        sys.exit(1)
    
    input_path = sys.argv[1]
    double_fasta_entry(input_path)


