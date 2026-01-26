import sys
import os

def concatenate_fasta(input_file):
    # Construct file paths
    path = os.path.dirname(input_file) or "."
    out_file = os.path.basename(input_file)
    output_fasta = os.path.join(path, f"concat.{out_file}")
    output_info = os.path.join(path, f"header_info.{out_file}.txt")
    header_prefix = f"{out_file}_flanks_numt"
    
    with open(input_file, 'r') as infile, open(output_fasta, 'w') as outfile_fasta, open(output_info, 'w') as outfile_info:
        header = ""
        sequence = ""
        count = 0
        id = 1
        
        for line in infile:
            if line.startswith(">"):
                count += 1
                if count % 2 == 1:
                    # Write the previous header and sequence
                    if header:
                        outfile_fasta.write(f">{header_prefix}{id}\n{sequence}\n")
                        outfile_info.write(f"{header_prefix}{id}\t{header}\n")
                        id += 1
                    header = line.strip()  # Start a new header
                    sequence = ""
                else:
                    header += line.strip()  # Append to the current header
            else:
                sequence += line.strip()  # Append to the sequence
        
        # Write the last header and sequence
        if header:
            outfile_fasta.write(f">{header_prefix}{id}\n{sequence}\n")
            outfile_info.write(f"{header_prefix}{id}\t{header}\n")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python concat_fasta.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    concatenate_fasta(input_file)


