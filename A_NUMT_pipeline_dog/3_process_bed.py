#!/usr/bin/env python3
import sys
import pandas as pd

# Minimal post-processing to normalize fields and add lengthAfterMerging

def process_bed_file(input_stream, output_stream):
    df = pd.read_csv(input_stream, sep='\t', header=None, names=['chrom', 'start', 'end', 'info', 'score', 'strand'])
    df['score'] = '.'

    def process_info(info):
        try:
            entries = info.split(',')
            max_length_entry = max(entries, key=lambda x: int(x.split('|')[1].split('=')[1]))
            processed_entries = [entry.replace('MT', 'MT_merged') if entry == max_length_entry else entry for entry in entries]
            return ','.join(processed_entries)
        except Exception:
            return info

    df['info'] = df['info'].apply(process_info)
    df['lengthAfterMerging'] = df.apply(lambda row: int(row['end']) - int(row['start']), axis=1)
    df['info'] = df.apply(lambda row: f"{row['info']}|lengthAfterMerging={row['lengthAfterMerging']}", axis=1)
    df.to_csv(output_stream, sep='\t', header=False, index=False)

if __name__ == '__main__':
    process_bed_file(sys.stdin, sys.stdout)
