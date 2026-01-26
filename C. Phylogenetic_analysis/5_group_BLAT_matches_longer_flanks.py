import pandas as pd
import os
import sys

list_ass = ['mGorGor1', 'mPanPan1', 'mPanTro3', 'mPonAbe1', 'mPonPyg2', 'mSymSyn1', 'CHM13']

def get_numt_loci(FLANK_SIZE):
    """
    Reads BED files for different query species, processes them to create a DataFrame
    with loci information for each NUMT, and returns the concatenated DataFrame.
    """
    list_df = []  # Initialize a list to hold DataFrames for each query species.
    
    # Iterate over the list of query species.
    for QUERY in list_ass:
        # Construct the file path for the current query species.
        numtDF = f"../numts/merged.blast.rawCoords.{QUERY}.bed"
        # Read the BED file into a DataFrame.
        df = pd.read_table(numtDF, header=None, names=['Chr', 'Start', 'End', 'Name', 'Score', 'Strand'])
        # Create a 'Label' column by combining Strand, Chr, Start, and End.
        df['Label'] = df['Strand'] + "::" + df['Chr'] + ":" + df['Start'].astype('str') + "-" + df['End'].astype('str')
        # Create a DataFrame with the relevant columns.
        iter_df = pd.DataFrame({
            'Query_Species': QUERY,
            'NUMT': QUERY + "_numt" + df['Score'].astype('str'),
            'NUMT_info': df['Name'],
            'Loci': df['Label'],
        })
        # Append the DataFrame to the list.
        list_df.append(iter_df)
    
    # Concatenate all DataFrames in the list into a single DataFrame.
    numt_loci = pd.concat(list_df)
    
    # Save the result to a tab-delimited file.
    numt_loci.to_csv(f'results_{FLANK_SIZE}bp/loci_numts.flanks{FLANK_SIZE}.tab', sep='\t', index=False)
    
    return numt_loci

def get_blat_matches(FLANK_SIZE):
    """
    Reads BED files for BLAT matches between different query and database species,
    processes them to create a DataFrame with match information, and returns the concatenated DataFrame.
    """
    list_df = []  # Initialize a list to hold DataFrames for each query and database species combination.

    # Iterate over the list of query species.
    for Q in list_ass:
        # Iterate over the list of database species.
        for DB in list_ass:
            # Construct the file path for the current query and database species combination.
            proc_blat_file = f"results_{FLANK_SIZE}bp/proc.{Q}.to.{DB}.flanks{FLANK_SIZE}.bed"

            # Check if the file exists
            if not os.path.exists(proc_blat_file):
                print(f"File not found: {proc_blat_file}")
                continue

            # Read the BED file into a DataFrame.
            df = pd.read_table(proc_blat_file, header=None, names=['Chr', 'Start', 'End', 'Name', 'Score', 'Strand'])

            # Check if the DataFrame is empty
            if df.empty:
                print(f"No data in file: {proc_blat_file}")
                continue

            # Create a 'Label' column by combining Strand, Chr, Start, and End.
            df['Label'] = df['Strand'] + "::" + df['Chr'] + ":" + df['Start'].astype('str') + "-" + df['End'].astype('str')

            # Group by 'Name' and aggregate unique 'Label' values as a list
            grouped_df = df.groupby("Name").agg({"Label": lambda x: sorted(x.unique())}).reset_index()

            # Check if grouped_df is empty
            if grouped_df.empty:
                print(f"No grouped data for {proc_blat_file}")
                continue

            # Remove the '_flanks1000' portion from the 'Name' column to match numt_loci format
            grouped_df['Name'] = grouped_df['Name'].str.replace(r'_flanks\d+', '', regex=True)

            # Create a DataFrame with the relevant columns.
            iter_df = pd.DataFrame({
                'Query_Species': Q,
                'NUMT': grouped_df['Name'],
                'DB_Species': DB,
                'Matches': grouped_df['Label']
            })

            # Append the DataFrame to the list.
            list_df.append(iter_df)

    # Check if list_df is empty before concatenation
    if not list_df:
        print("No data to concatenate.")
        return pd.DataFrame()

    # Concatenate all DataFrames in the list into a single DataFrame.
    blat_matches = pd.concat(list_df).reset_index(drop=True)

    # Get NUMT loci information.
    numt_loci = get_numt_loci(FLANK_SIZE=FLANK_SIZE)

    # Check if numt_loci is empty
    if numt_loci.empty:
        print("NUMT loci data is empty.")
        return blat_matches

    # Merge the NUMT loci information with BLAT matches.
    blat_matches = pd.merge(numt_loci, blat_matches)

    # Count the number of matches per entry.
    blat_matches['Counts'] = blat_matches['Matches'].apply(len)

    # Check if blat_matches is empty before saving
    if blat_matches.empty:
        print("Final DataFrame is empty, nothing to save.")
    else:
        # Save the result to a tab-delimited file.
        blat_matches.to_csv(f"results_{FLANK_SIZE}bp/blat_matches.flanks{FLANK_SIZE}.tab", sep='\t', index=False)

    return blat_matches


if __name__ == "__main__":
    # Ensure the user has provided a FLANK_SIZE argument
    if len(sys.argv) != 2:
        print("Usage: python script.py <FLANK_SIZE>")
        sys.exit(1)

    # Get the FLANK_SIZE from the command-line arguments
    FLANK_SIZE = int(sys.argv[1])

    # Run the function and store the result in a DataFrame.
    get_blat_matches(FLANK_SIZE)
