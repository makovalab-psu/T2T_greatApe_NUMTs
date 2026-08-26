import sys
import pandas as pd
import numpy as np


def blast_to_bed(input_file):
    # Define the column names from the BLAST output, replacing spaces with underscores
    column_names = [
        "query_acc.ver", "subject_acc.ver", "percent_identity", "alignment_length",
        "mismatches", "gap_opens", "q_start", "q_end", "s_start",
        "s_end", "evalue", "bit_score"
    ]
    
    # Read the file, filtering out lines that start with '#'
    with open(input_file, 'r') as file:
        data_lines = [line.strip() for line in file if not line.startswith('#')]
    
    # Create a DataFrame from the filtered data
    try:
        df = pd.DataFrame([line.split('\t') for line in data_lines], columns=column_names)
    except ValueError as e:
        print(f"Error creating DataFrame: {e}")
        sys.exit(1)
    
    # Determine strand based on s_start and s_end
    def determine_strand(row):
        return '+' if int(row['s_start']) < int(row['s_end']) else '-'
    
    # Function to flip start and end coordinates if on the negative strand
    def adjust_coordinates(row):
        if row['strand'] == '-':
            return pd.Series({
                'chromStart': int(row['s_end']) - 1,  # Convert to 0-based
                'chromEnd': int(row['s_start'])
            })
        else:
            return pd.Series({
                'chromStart': int(row['s_start']) - 1,  # Convert to 0-based
                'chromEnd': int(row['s_end'])
            })
    
    # Convert to BED format (subject, start, end, name, score, strand)
    bed_df = pd.DataFrame()
    bed_df['chrom'] = df['subject_acc.ver']
    df['strand'] = df.apply(determine_strand, axis=1)
    bed_df[['chromStart', 'chromEnd']] = df.apply(adjust_coordinates, axis=1)

    # Reduce coordinates that are above the original breakpoint (mtDNA size) in the mtDNA doubled genome.
    mt_size = list_mtdna[assembly]
    # Ensure 'q_start' and 'q_end' are of type integer and adjust based on the mtDNA size.
    df['q_start'] = df['q_start'].astype('int')
    df['q_start'] = np.where(df['q_start'] > mt_size, df['q_start'] - mt_size, df['q_start'])
    df['q_end'] = df['q_end'].astype('int')
    df['q_end'] = np.where(df['q_end'] > mt_size, df['q_end'] - mt_size, df['q_end'])

    # Adding fragment length and percent identity to the name column
    bed_df['name'] = df.apply(
        lambda row: f"{row['query_acc.ver']}:{row['q_start']}-{row['q_end']}|length={bed_df.loc[row.name, 'chromEnd'] - bed_df.loc[row.name, 'chromStart']}|identity={row['percent_identity']}|nuclearStrand={row['strand']}",
        axis=1
    )
    
    bed_df['score'] = '.'  # Convert score to integer
    bed_df['strand'] = df['strand']
    
    # Output to stdout
    bed_df.to_csv(sys.stdout, sep='\t', header=False, index=False)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python blast_to_bed.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # MtDNA sizes per species.
    list_mtdna = { 'hg002':16569, 'mGorGor1':16407, 'chm13':16569, 'CHM13':16569, 'mPanPan1':16569, 'mPanTro3':16620, 'mPonAbe1':16496, 'mPonPyg2':16461, 'mSymSyn1':16515 }

    assembly = None
    for x in list_mtdna.keys():
        # print(x,input_file)
        if x in input_file:
            assembly = x
            break
        else:
            continue
    if not assembly:
        sys.exit("# Cannot find assembly name in the input file name.")
    
    blast_to_bed(input_file)
