import sys
import pandas as pd

def filter_by_dissimilarity(kmer_dissim_threshold,
                            input_file_path='expanded_matches_similarity_results.tab',
                            output_file_path='kmer_test/kmerFilter_expanded_matches_similarity_results.tab'):
    """
    Reads a tab-delimited input file, filters the DataFrame to include only rows where 'Hit_Dissimilarity' is greater than 
    or equal to the specified threshold, and saves the filtered DataFrame to a tab-delimited file.

    Parameters:
    - input_file_path (str): The path to the input tab-delimited file. Default is 'expanded_matches_similarity_results.tab'.
    - threshold (float): The threshold value for 'Hit_Dissimilarity'. Default is 0.5.
    - output_file_path (str): The file path to save the filtered DataFrame. Default is 'filtered_data.tab'.

    Returns:
    - DataFrame: A filtered DataFrame containing only rows with 'Hit_Dissimilarity' >= threshold.
    """
    # Read the input data file.
    try:
        df = pd.read_table(input_file_path)
    except:
        print(f"### File or path does not exist: {input_file_path}")
        raise ValueError
    
    # Filter the DataFrame based on the 'Hit_Dissimilarity' column.
    df_filtered = df[df['Hit_Dissimilarity'] <= kmer_dissim_threshold]
    
    # Save the filtered DataFrame to a tab-delimited file.
    df_filtered.to_csv(output_file_path, sep='\t', index=False)
    
    return df_filtered


def create_species_matrix(  kmer_dissim_threshold, 
                            input_file_path, 
                            output_file_path ):
    """
    Transforms the input DataFrame into a matrix where each row represents a NUMT,
    and each column represents a species. The matrix indicates the presence (1) or 
    absence (0) of species for each NUMT, and saves the matrix to a file.

    Parameters:
    - df_sim (DataFrame): Input DataFrame with columns 'NUMT' and 'DB_CommonName'.
    - output_file_path (str): Path to save the resulting matrix file. Default is 'species_matrix.tsv'.

    Returns:
    - DataFrame: A matrix where rows are NUMTs and columns are species.
    """
    # Make a copy of the input DataFrame
    df = filter_by_dissimilarity(kmer_dissim_threshold, input_file_path=input_file_path, output_file_path=output_file_path )

    # Remove Siamang.
    df = df[(df['DB_CommonName']!='Siamang')&(df['Query_CommonName']!='Siamang')]
    
    # Create a list to hold the processed data
    processed_data = []

    # Iterate over each unique NUMT
    for numt in df['NUMT'].drop_duplicates():
        # Get the unique species names for the current NUMT
        species_names = df[df['NUMT'] == numt]['DB_CommonName'].drop_duplicates().tolist()
        # Append the data to the list
        processed_data.append([numt, species_names, len(species_names)])

    # Create a DataFrame from the processed data
    df_processed = pd.DataFrame(processed_data, columns=['NUMT', 'List_DB_CommonName', 'Counts'])

    # Extract a sorted list of unique species from the List_DB_CommonName column
    #unique_species = sorted(set(species for sublist in df_processed['List_DB_CommonName'] for species in sublist))

    unique_species = [ "Bonobo", "Chimpanzee", "Human", "Gorilla", "Sorang", "Borang" ]  #, "Siamang"]

    # Initialize a DataFrame to hold the species matrix
    species_matrix = pd.DataFrame(0, index=df_processed['NUMT'], columns=unique_species)

    # Populate the species matrix
    for idx, row in df_processed.iterrows():
        for species in row['List_DB_CommonName']:
            species_matrix.at[row['NUMT'], species] = 1

    # Save the resulting matrix to a file
    species_matrix.to_csv(output_file_path, sep='\t', index=True)
    print(f"Matrix saved to {output_file_path}")

    # Write the Nexus file
    nexus_file_name = output_file_path.replace('.tsv', '.nexus')
    write_nexus_file(species_matrix, output_file_path=nexus_file_name)

    # Display the resulting matrix
    return species_matrix


def write_nexus_file(matrix, output_file_path):
    nexus_file = open(output_file_path, "w") # open the nexus file

    n_taxa = 6 #excludes Siamang
    # n_taxa = 7 #includes Siamang

    # taxonomy in nexus file
    nexus_file.write("#NEXUS\n\n") 
    # nexus_file.write("[begin taxa;]\n")
    # nexus_file.write(f"[dimensions ntax={n_taxa} nchar={len(matrix)};]\n")
    # nexus_file.write("[taxlabels]\n")
    # nexus_file.write("[Bonobo]\n")
    # nexus_file.write("[Chimpanzee]\n")
    # nexus_file.write("[Human]\n")
    # nexus_file.write("[Gorilla]\n")
    # nexus_file.write("[B._orangutan]\n")
    # nexus_file.write("[S._orangutan]\n")
    # # nexus_file.write("[Siamang]\n")
    # nexus_file.write("[end;]\n\n")

    # data of the nexus file
    nexus_file.write(f"Begin data;\n")
    nexus_file.write(f"dimensions ntax={n_taxa} nchar={len(matrix)};\n")
    nexus_file.write(f"matrix\n")
    nexus_file.write(f"Bonobo {''.join(list(str(int(i)) for i in matrix['Bonobo']))}\n")
    nexus_file.write(f"Chimpanzee {''.join(list(str(int(i)) for i in matrix['Chimpanzee']))}\n")
    nexus_file.write(f"Human {''.join(list(str(int(i)) for i in matrix['Human']))}\n")
    nexus_file.write(f"Gorilla {''.join(list(str(int(i)) for i in matrix['Gorilla']))}\n")
    nexus_file.write(f"B._orangutan {''.join(list(str(int(i)) for i in matrix['Borang']))}\n")
    nexus_file.write(f"S._orangutan {''.join(list(str(int(i)) for i in matrix['Sorang']))}\n")
    # nexus_file.write(f"Siamang {''.join(list(str(int(i)) for i in matrix['Siamang']))}\n")
    nexus_file.write(";\n")
    nexus_file.write("end;\n")

    nexus_file.close() #close the nexus file
    print(f"Matrix saved to {output_file_path}")


if __name__ == "__main__":
    # Ensure the correct number of arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <FLANK_SIZE>")
        sys.exit(1)

    # Parse the kmer_dissim_threshold argument
    try:
        FLANK_SIZE = int(sys.argv[1])
        list_kd = [ 0.2, 0.5, 0.8 ]
        #FLANK_SIZE=600
        DIR=f"results_{FLANK_SIZE}bp/"
        for kmer_dissim_threshold in list_kd:
            create_species_matrix(  kmer_dissim_threshold=kmer_dissim_threshold, 
                                    input_file_path=f'{DIR}expanded_matches_similarity_results.flanks{FLANK_SIZE}.tab',
                                    output_file_path=f'{DIR}species_matrix_flanks{FLANK_SIZE}_{kmer_dissim_threshold}.tsv' )
    except:
        pass
