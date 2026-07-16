import sys
import pandas as pd

def get_matrix(sample, input_file_path):
    input_file_path = f'combined_matrix.t2t_{sample}.tsv'
    df = pd.read_table( input_file_path )
    # Values that can be inferred form the rest of the table.
    df = df.fillna('?')

    # For rows where all T2T columns are 'Missing', set the matching species to 1, and the non-matching species columns to 0.
    species_cols = ['Bonobo', 'Chimpanzee', 'Human', 'Gorilla', 'Sorang', 'Borang']
    all_missing_mask = df[species_cols].eq('?').all(axis=1)
    numt_to_species = {
		'CHM13': 'Human',
		'mGorGor': 'Gorilla',
		'mPanPan': 'Bonobo',
		'mPanTro': 'Chimpanzee',
		'mPonPyg': 'Sorang',
		'mPonAbe': 'Borang'
	}
    for idx in df[all_missing_mask].index:
       numt_id = df.loc[idx, 'NUMT_ID']
       matching_species = next((sp for prefix, sp in numt_to_species.items() if numt_id.startswith(prefix)), None)
       if matching_species:
           df.loc[idx, species_cols] = 0
           df.loc[idx, matching_species] = 1

    # Drop rows if they still are missing ('?').
    df = df[~df.isin(['?']).any(axis=1)]
    return df


def write_nexus_file(matrix, output_file_path):
    nexus_file = open(output_file_path, "w") # open the nexus file

    n_taxa = len(matrix.columns) - 1
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
        if col not in ['NUMT_ID','Bonobo', 'Chimpanzee', 'Human', 'Gorilla', 'Sorang', 'Borang']:
            nexus_file.write(f"{col} {''.join(list(str(int(i)) for i in matrix[col]))}\n")
        else:
            pass
    nexus_file.write(";\n")
    nexus_file.write("end;\n")

    nexus_file.close() #close the nexus file
    print(f"Matrix saved to {output_file_path}")

if __name__ == "__main__":
    # # Ensure the correct number of arguments is provided
    # if len(sys.argv) != 2:
    #     print("Usage: python script.py <FLANK_SIZE>")
    #     sys.exit(1)

    # # Parse the kmer_dissim_threshold argument
    # try:
        # SAMPLE = int(sys.argv[1])
        list_samples = [ 'hprc', 'GorGor', 'PanPan', 'PanTro', 'PonAbe', 'PonPyg' ]
        for SAMPLE in list_samples:
            matrix = get_matrix( sample=SAMPLE, input_file_path=f'combined_matrix.t2t_{SAMPLE}.tsv' )
            write_nexus_file( matrix, output_file_path=f'combined_matrix.t2t_{SAMPLE}.nexus' )
    # except:
    #     print('# What happened?')
    #     pass
