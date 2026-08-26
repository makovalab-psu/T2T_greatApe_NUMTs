#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
import os

# This script converts BLAST -outfmt 7 output to BED.
# It uses mt_sizes.tsv (tab: assembly\tmt_size) in the pipeline folder to get mt length.

def read_mt_sizes(mt_sizes_file):
    sizes = {}
    if not os.path.exists(mt_sizes_file):
        return sizes
    with open(mt_sizes_file) as fh:
        for line in fh:
            line=line.strip()
            if not line or line.startswith('#'): continue
            parts=line.split()  # tab or space
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def blast_to_bed(input_file, mt_sizes):
    # Read relevant non-comment lines
    data_lines = []
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip()=='' or line.startswith('Query='):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 12:
                data_lines.append(parts[:12])

    if not data_lines:
        return pd.DataFrame(columns=['chrom','chromStart','chromEnd','name','score','strand'])

    column_names = [
        "query_acc.ver", "subject_acc.ver", "percent_identity", "alignment_length",
        "mismatches", "gap_opens", "q_start", "q_end", "s_start",
        "s_end", "evalue", "bit_score"
    ]
    df = pd.DataFrame(data_lines, columns=column_names)

    # Determine assembly by searching known keys in input filename
    assembly = None
    for k in mt_sizes.keys():
        if k in os.path.basename(input_file):
            assembly = k
            break
    # fallback: try any assembly snippet in file name
    if not assembly:
        # pick first available size if only one specified
        if len(mt_sizes) == 1:
            assembly = list(mt_sizes.keys())[0]

    # Determine strand
    df['strand'] = df.apply(lambda r: '+' if int(r['s_start']) < int(r['s_end']) else '-', axis=1)

    def adjust_coords(row):
        if row['strand'] == '-':
            start = int(row['s_end']) - 1
            end = int(row['s_start'])
        else:
            start = int(row['s_start']) - 1
            end = int(row['s_end'])
        return pd.Series({'chromStart': start, 'chromEnd': end})

    bed_df = pd.DataFrame()
    bed_df['chrom'] = df['subject_acc.ver']
    bed_df[['chromStart','chromEnd']] = df.apply(adjust_coords, axis=1)

    # adjust q coordinates using mt size if available
    mt_size = mt_sizes.get(assembly, None)
    df['q_start'] = df['q_start'].astype(int)
    df['q_end'] = df['q_end'].astype(int)
    if mt_size:
        df['q_start'] = np.where(df['q_start'] > mt_size, df['q_start'] - mt_size, df['q_start'])
        df['q_end'] = np.where(df['q_end'] > mt_size, df['q_end'] - mt_size, df['q_end'])

    bed_df['name'] = df.apply(lambda row: f"{row['query_acc.ver']}:{row['q_start']}-{row['q_end']}|length={int(row['alignment_length'])}|identity={row['percent_identity']}|nuclearStrand={row['strand']}", axis=1)
    bed_df['score'] = '.'
    bed_df['strand'] = df['strand']

    return bed_df[['chrom','chromStart','chromEnd','name','score','strand']]


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: 2_process_blast.py <blast_outfmt7_file>', file=sys.stderr)
        sys.exit(1)
    infile = sys.argv[1]
    mt_sizes = read_mt_sizes(os.path.join(os.path.dirname(__file__), 'mt_sizes.tsv'))
    bed = blast_to_bed(infile, mt_sizes)
    bed.to_csv(sys.stdout, sep='\t', header=False, index=False)
