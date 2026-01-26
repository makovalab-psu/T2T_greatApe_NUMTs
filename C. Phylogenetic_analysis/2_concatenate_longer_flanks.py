import sys

def concatenate_fasta(assembly,flank_size=None):
    if not assembly:
        raise ValueError("Assembly argument is required.")
    
    input_file = f"flanks{flank_size}.{assembly}.fasta" if flank_size else f"flanks.{assembly}.fasta"
    output_fasta = f"concat.flanks{flank_size}.{assembly}.fasta" if flank_size else f"concat.flanks.{assembly}.fasta"
    output_info = f"header_info.{assembly}.flanks{flank_size}.txt"
    header_prefix = f"{assembly}_flanks{flank_size}_numt" if flank_size else f"{assembly}_flanks_numt"
    
    with open(input_file, 'r') as infile, open(output_fasta, 'w') as outfile_fasta, open(output_info, 'w') as outfile_info:
        header = ""
        sequence = ""
        count = 0
        id = 1
        
        for line in infile:
            if line.startswith(">"):
                count += 1
                if count % 2 == 1:
                    if header:
                        outfile_fasta.write(f">{header_prefix}{id}\n{sequence}\n")
                        outfile_info.write(f"{header_prefix}{id}\t{header}\n")
                        id += 1
                    header = line.strip()
                    sequence = ""
                else:
                    header += line.strip()
            else:
                sequence += line.strip()
        
        if header:
            outfile_fasta.write(f">{header_prefix}{id}\n{sequence}\n")
            outfile_info.write(f"{header_prefix}{id}\t{header}\n")

if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     raise ValueError("Assembly argument is required.")
    # assembly = sys.argv[1]

    flank_size='500'
    for assembly in ['mGorGor1', 'mPanPan1', 'mPanTro3', 'mPonAbe1', 'mPonPyg2', 'mSymSyn1', 'CHM13']:
        concatenate_fasta(assembly,flank_size)
