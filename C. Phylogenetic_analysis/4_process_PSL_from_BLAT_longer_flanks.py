import pandas as pd
import sys

def process_psl(FLANK_SIZE):
    header_psl = ['match', 'mismatch', 'rep_name', 'Ns', 'Q_gap_counts', 'Q_gap_bases', 'T_gap_counts', 'T_gap_bases', 'strand', 'Q_name', 
                  'Q_size', 'Q_start', 'Q_end', 'T_name', 'T_size', 'T_start', 'T_end']

    # Only the first 17 columns for the BLAT PSL are used by Mark Loftus to consider the quality of DB hit.
    df = pd.read_table(INPUT_FILE, skiprows=5, header=None, names=header_psl, usecols=range(17))

    # Handle strand information for query and target
    df['Q_strand'] = df['strand'].apply(lambda x: x[0])
    df['T_strand'] = df['strand'].apply(lambda x: x[-1] if len(x) > 1 else '+')

    # Annotate with Numt coordinates
    df_numt = pd.read_table(f"results_{FLANK_SIZE}bp/header_info.{QUERY}.flanks{FLANK_SIZE}.txt", header=None, names=['Q_name', 'Info'])
    df = df.merge(df_numt)

    # NUMT coordinates
    df['Numt_chr'] = df['Info'].str.split('>').str[-2].str.split(':').str[2]
    df['Numt_start'] = df['Info'].str.split('>').str[-2].str.split(':').str[3].str.split('-').str[1].astype('int')
    df['Numt_end'] = df['Info'].str.split('>').str[-1].str.split(':').str[3].str.split('-').str[0].astype('int')

    # Calculate Numt size
    df['Numt_size'] = df['Info'].str.split('_').str[1].str[:-2].astype('int')

    # Alternate Numt size
    df['Numt_size2'] = df['Numt_end'] - df['Numt_start']
    negative_size_mask = df['Numt_size2'] < 0
    df.loc[negative_size_mask, ['Numt_start', 'Numt_end']] = df.loc[negative_size_mask, ['Numt_end', 'Numt_start']].values
    df['Numt_size2'] = df['Numt_end'] - df['Numt_start']

    # Keep matches with a total gap size (between flanks) within +/- 30% of the NUMT size
    df["Diff_size"] = abs(df['Numt_size'] - df['T_gap_bases']) / df['Numt_size']
    df_within30 = df[df["Diff_size"] < 0.3]

    # Create BED columns from PSL columns
    df_bed = pd.DataFrame({
        "T_name": df_within30["T_name"],
        "T_start": df_within30["T_start"],
        "T_end": df_within30["T_end"],
        "Q_name": df_within30["Q_name"],
        "Empty_Score": 0,  # Default score value
        "T_strand": df_within30["T_strand"]
    })

    # Export BED files of DB hits
    df_bed.to_csv(OUT_BED, sep='\t', index=None, header=None)

if __name__ == "__main__":
    # Ensure the user has provided a FLANK_SIZE argument
    if len(sys.argv) != 2:
        print("Usage: python 4_process_PSL_from_BLAT_longer_flanks.py <FLANK_SIZE>")
        sys.exit(1)

    # Get the FLANK_SIZE from the command-line arguments
    FLANK_SIZE = int(sys.argv[1])
    DIR=f"results_{FLANK_SIZE}bp/"
    for QUERY in ['mGorGor1', 'mPanPan1', 'mPanTro3', 'mPonAbe1', 'mPonPyg2', 'mSymSyn1', 'CHM13']:
        for DB in ['mGorGor1', 'mPanPan1', 'mPanTro3', 'mPonAbe1', 'mPonPyg2', 'mSymSyn1', 'CHM13']:
            INPUT_FILE = f"{DIR}out.{QUERY}.to.{DB}.flanks{FLANK_SIZE}.psl"
            OUT_FILE = f"{DIR}proc.{QUERY}.to.{DB}.flanks{FLANK_SIZE}.tab"
            OUT_BED = f"{DIR}proc.{QUERY}.to.{DB}.flanks{FLANK_SIZE}.bed"
            # Run the function
            process_psl(FLANK_SIZE)
