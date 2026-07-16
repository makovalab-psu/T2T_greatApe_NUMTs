import sys
import pandas as pd

ERROR_FILE_NOT_FOUND = 1
ERROR_PROCESSING = 2

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
    except FileNotFoundError:
        print(f"ERROR {ERROR_FILE_NOT_FOUND}: file not found: {input_file_path}")
        raise
    except Exception as exc:
        print(f"ERROR {ERROR_PROCESSING}: failed to read input file {input_file_path}: {exc}")
        raise
    
    # Filter the DataFrame based on the 'Hit_Dissimilarity' column.
    df_filtered = df[df['Hit_Dissimilarity'] <= kmer_dissim_threshold]
    
    # Save the filtered DataFrame to a tab-delimited file.
    df_filtered.to_csv(output_file_path, sep='\t', index=False)
    
    # print(df_filtered)
    return df_filtered


def create_species_matrix(  kmer_dissim_threshold, 
                            input_file_path, 
                            t2t_file,
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
    # print("Filtered df size:", df.shape[0])

    # Remove Siamang.
    df['DB_CommonName'] = df['DB_Species']
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
    unique_species = sorted(set(species for sublist in df_processed['List_DB_CommonName'] for species in sublist))
    # print(unique_species)

    # unique_species = [ "Bonobo", "Chimpanzee", "Human", "Gorilla", "Sorang", "Borang" ]  #, "Siamang"]

    # Initialize a DataFrame to hold the species matrix
    species_matrix = pd.DataFrame(0, index=df_processed['NUMT'], columns=unique_species)

    # Populate the species matrix
    for idx, row in df_processed.iterrows():
        for species in row['List_DB_CommonName']:
            species_matrix.at[row['NUMT'], species] = 1

    # Add P/A matrix results for the T2T asssemblies.
    species_matrix = add_t2t_matrix( species_matrix, t2t_file )

    # Save the resulting matrix to a file
    species_matrix.to_csv(output_file_path, sep='\t', index=True)
    print(f"Matrix saved to {output_file_path}")

    # Write the Nexus file
    nexus_file_name = output_file_path.replace('.tsv', '.nexus')
    write_nexus_file(species_matrix, output_file_path=nexus_file_name)

    # Display the resulting matrix
    return species_matrix

def add_t2t_matrix(matrix, t2t_file):
    t2t = pd.read_table(t2t_file).set_index('NUMT')
    # print(t2t)
    matrix2 = pd.merge(t2t, matrix, left_index=True, right_index=True, how='outer').fillna(0)
    matrix2 = matrix2.astype(int)
    print(matrix2)
    return matrix2


def write_nexus_file(matrix, output_file_path):
    nexus_file = open(output_file_path, "w") # open the nexus file

    n_taxa = len(matrix.columns)
    print(n_taxa)
    # n_taxa = 6 #excludes Siamang
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
    for col in matrix.columns:
        if col not in ['NUMT','Bonobo', 'Chimpanzee', 'Human', 'Gorilla', 'Sorang', 'Borang']:
            nexus_file.write(f"{col} {''.join(list(str(int(i)) for i in matrix[col]))}\n")
        else:
            pass
    nexus_file.write(";\n")
    nexus_file.write("end;\n")

    nexus_file.close() #close the nexus file
    print(f"Matrix saved to {output_file_path}")


if __name__ == "__main__":
    # Parse the kmer_dissim_threshold argument
    try:
        FLANK_SIZE = 500
        kmer_dissim_threshold = 0.5
        create_species_matrix(  kmer_dissim_threshold=kmer_dissim_threshold, 
                                input_file_path=f'expanded_matches_similarity_results.flanks{FLANK_SIZE}.tab',
                                t2t_file='t2t.species_matrix_flanks500_0.5.tsv',
                                output_file_path=f'hprc5.species_matrix_flanks{FLANK_SIZE}_{kmer_dissim_threshold}.tsv' )
    except FileNotFoundError:
        sys.exit(ERROR_FILE_NOT_FOUND)
    except Exception as exc:
        print(f"ERROR {ERROR_PROCESSING}: {exc}")
        sys.exit(ERROR_PROCESSING)
