import sys
import pandas as pd

def process_bed_file(input_stream, output_stream):
    # Read the BED file from stdin
    df = pd.read_csv(input_stream, sep='\t', header=None, names=['chrom', 'start', 'end', 'info', 'score', 'strand'])
    
    # Replace 'score' column with periods
    df['score'] = '.'
    
    # Function to extract the largest length value and replace 'MT' with 'MT_merged'
    def process_info(info):
        try:
            entries = info.split(',')
            # Extract length from each entry and find the maximum
            max_length_entry = max(entries, key=lambda x: int(x.split('|')[1].split('=')[1]))
            processed_entries = [entry.replace('MT', 'MT_merged') if entry == max_length_entry else entry for entry in entries]
            return ','.join(processed_entries)
        except IndexError:
            sys.exit()

    # Apply the function to each row
    df['info'] = df['info'].apply(process_info)
    
    # Add lengthAfterMerging to the info column
    df['lengthAfterMerging'] = df.apply(lambda row: row['end'] - row['start'], axis=1)
    df['info'] = df.apply(lambda row: row['info'] + f'|lengthAfterMerging={row["lengthAfterMerging"]}', axis=1)
    
    # Write the processed DataFrame to stdout
    df.to_csv(output_stream, sep='\t', header=False, index=False)

if __name__ == "__main__":
    process_bed_file(sys.stdin, sys.stdout)
