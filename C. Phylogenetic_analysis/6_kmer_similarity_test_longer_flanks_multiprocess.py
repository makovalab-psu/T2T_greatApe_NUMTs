import pandas as pd
import ast
import pysam
from Bio.Seq import Seq
from tqdm import tqdm
from multiprocessing import Pool, cpu_count


FLANK_SIZE = 750


# Function to process a single row
def process_row(row):
    mykmerSize = 14
    dissimilarities = []
    similarities = []

    # Extract the contig, start, and end positions from the 'Loci' column
    contig = str(row['Loci'].split(":")[2])
    start = str(int(row['Loci'].split(":")[3].split("-")[0]) - 500)
    end = str(int(row['Loci'].split(":")[3].split("-")[1]) + 500)

    # Determine the focus species
    focus_species = str(row['Query_Species'])

    try:
        # Fetch the reference sequence from the genome
        sequence = ''.join(pysam.faidx(GenomeDict2[focus_species], f"{contig}:{start}-{end}").split()[1:])
    except Exception as e:
        print(f"Error fetching reference sequence: {e}")
        return [], []  # Return empty lists if the sequence cannot be fetched

    # Process each match in the 'Matches' column
    for match in ast.literal_eval(str(row['Matches'])):
        querySpecies = str(row['DB_Species'])

        # Determine the coordinate for the query sequence
        coordinate = str(match.split("::")[1])

        try:
            # Fetch and possibly reverse complement the query sequence
            if str(match.split("::")[0]) == '+':
                sequence2 = ''.join(pysam.faidx(GenomeDict2[querySpecies], coordinate).split()[1:])
            else:
                presequence = Seq(''.join(pysam.faidx(GenomeDict2[querySpecies], coordinate).split()[1:]))
                sequence2 = str(presequence.reverse_complement())
        except Exception as e:
            print(f"Error fetching query sequence: {e}")
            continue  # Skip this iteration if the sequence cannot be fetched

        # Generate k-mers and insertion k-mers
        kmers = []
        insertionkmers = []
        i = 0
        while i < (len(sequence) - (mykmerSize - 1)):
            kmers.append(str(sequence[i:i + mykmerSize]).upper())
            insertionkmers.append(str(sequence2[i:i + mykmerSize]).upper())
            i += 1

        # Create a DataFrame to hold k-mer counts
        uniqueKmers = set(kmers + insertionkmers)
        tempDF = pd.DataFrame(0, index=['Reference', 'MEI'], columns=[x for x in uniqueKmers])

        # Count the occurrences of k-mers
        for kmer in kmers:
            tempDF.at['Reference', kmer] += 1
        for kmer2 in insertionkmers:
            tempDF.at['MEI', kmer2] += 1

        # Calculate the differences and dissimilarities
        tempDF.loc['Diff'] = tempDF.loc['Reference'] - tempDF.loc['MEI']
        tempDF.loc['Sum'] = tempDF.loc['Reference'] + tempDF.loc['MEI']
        dissimilarities.append(sum(abs(tempDF.loc['Diff'])) / sum(abs(tempDF.loc['Sum'])))
        similarities.append(1 - sum(abs(tempDF.loc['Diff'])) / sum(abs(tempDF.loc['Sum'])))

    return dissimilarities, similarities


# Read the tab-delimited file into a DataFrame
dataDF = pd.read_table(f"blat_matches.flanks{FLANK_SIZE}.tab")


# Dictionary mapping species to their corresponding reference genome FASTA files
GenomeDict2 = {
    "mGorGor1": "../../data/refs/mGorGor1.pri.cur.20231122.fasta",
    "mPanPan1": "../../data/refs/mPanPan1.pri.cur.20231122.fasta",
    "mPanTro3": "../../data/refs/mPanTro3.pri.cur.20231122.fasta",
    "mPonAbe1": "../../data/refs/mPonAbe1.pri.cur.20231205.fasta",
    "mPonPyg2": "../../data/refs/mPonPyg2.pri.cur.20231122.fasta",
    "mSymSyn1": "../../data/refs/mSymSyn1.pri.cur.20231205.fasta",
    "CHM13": "../../data/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
}

# Prepare the rows for processing
rows = [dataDF.loc[i] for i in dataDF.index]

# Use multiprocessing to process rows in parallel
with Pool(cpu_count()) as pool:
    results = list(tqdm(pool.imap(process_row, rows), total=len(rows), desc="Processing rows", unit="row"))

# Separate the results into two lists: seqDiss and seqSims
seqDiss, seqSims = zip(*results)

# Add the similarity results to the DataFrame and save to a CSV file
dataDF['Hit_Dissimilarity'] = seqDiss
dataDF['Hit_Similarity'] = seqSims
dataDF.to_csv(f"matches_similarity_results.flanks{FLANK_SIZE}.tab", index=None, sep='\t')


def expand_columns(file_path=f'matches_similarity_results.flanks{FLANK_SIZE}.tab',
                   output_file_path=f'expanded_matches_similarity_results.flanks{FLANK_SIZE}.tab'):
    """
    Expands rows of the DataFrame based on the lists in 'Matches', Hit_Dissimilarity and 'Hit_Similarity' columns.
    
    Parameters:
    - file_path (str): Path to the input tab-delimited file.
    - output_file_path (str): Path to save the expanded DataFrame.

    Returns:
    - DataFrame: The expanded DataFrame.
    """
    # Load the data from the file
    df = pd.read_table(file_path)

    # Add common species name.
    df['Query_CommonName'] = [{
        'CHM13':'Human', 'mGorGor1':'Gorilla', 'mPanPan1':'Bonobo', 'mPanTro3':'Chimpanzee', 'mSymSyn1':'Siamang', 'mPonAbe1':'Sorang', 'mPonPyg2':'Borang'
    }[x] for x in df['Query_Species'] ]
    df['DB_CommonName'] = [{
        'CHM13':'Human', 'mGorGor1':'Gorilla', 'mPanPan1':'Bonobo', 'mPanTro3':'Chimpanzee', 'mSymSyn1':'Siamang', 'mPonAbe1':'Sorang', 'mPonPyg2':'Borang'
    }[x] for x in df['DB_Species'] ]
    
    # Convert 'Matches' and 'Hit_Dissimilarity' from string representation of list to actual lists
    df['Matches'] = df['Matches'].apply(ast.literal_eval)
    df['Hit_Dissimilarity'] = df['Hit_Dissimilarity'].apply(ast.literal_eval)
    df['Hit_Similarity'] = df['Hit_Similarity'].apply(ast.literal_eval)

    rows = []
    for i, row in df.iterrows():
        matches = row['Matches']
        dissim = row['Hit_Dissimilarity']
        sim = row['Hit_Similarity']
        if len(matches) != len(dissim):
            print(matches[0])
            print(len(matches), len(dissim), len(sim))
            raise Exception("### The number of values do not match between Matches and Hit_Dissimilarity ###")
        else:
            for j in range(len(matches)):
                new_row = row.copy()
                new_row['Matches'] = matches[j]
                new_row['Hit_Dissimilarity'] = dissim[j]
                new_row['Hit_Similarity'] = sim[j]
                new_row['List_Index'] = j
                rows.append(new_row)
    df_expanded = pd.DataFrame(rows)
    
    # Save the expanded DataFrame to a new file
    df_expanded.to_csv(output_file_path, sep='\t', index=False)
    return df_expanded


# Call the function and get the expanded DataFrame
expand_columns(file_path=f'matches_similarity_results.flanks{FLANK_SIZE}.tab',
               output_file_path=f'expanded_matches_similarity_results.flanks{FLANK_SIZE}.tab')

