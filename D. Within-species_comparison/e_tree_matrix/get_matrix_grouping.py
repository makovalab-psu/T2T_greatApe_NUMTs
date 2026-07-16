import pandas as pd
import glob
import sys

dict_samples = { 'PanPan':'Bonobo', 'PanTro':'Chimpanzee', 'hprc':'Human', 
				 'GorGor':'Gorilla', 'PonAbe':'Sorang', 'PonPyg':'Borang' }

def get_matrix(sample):
	df = pd.read_table( f'combined_matrix.t2t_{sample}.tsv', dtype=str )
	df = df.replace({"1.0": "1", "0.0": "0"}) #make sure these are not saved as floats
	return df

def classify_values(s):
	if s != s:  # NaN check
		return 'NaN'
	ones = s.count('1')
	zeros = s.count('0')
	total = ones + zeros
	# Classifications.
	if zeros == 0:
		return 'Fixed' 	# Only 1s
	elif ones == 0:
		return 'Only in T2T' 	# Only 0s
	elif zeros == 1:
		return 'Almost fixed' 	# Mostly 1s and one 0
	elif ones == 1:
		return 'Rare'	# Mostly 0s and one 1
	else:
		return 'Variable' 	# Mix of 1s and 0s
	
def get_species_specific(sample, quiet=False):
	# Get matrix, then species and assembly names.
	df = get_matrix(sample)
	species = dict_samples[sample]
	assembly = df['NUMT_ID'].str.split('_').str[0][0]
	if not quiet:
		print( f"# Number of '{assembly}' NUMTs:", df.shape[0])
		
	# Get the order of T2T columns in this matrix.
	t2t_cols = list(df.columns[1:7])
	if not quiet:
		print( f"# The order of T2T columns is:", t2t_cols)
	
	# Unit test: are these the correct and complete valeus for the T2T columns?
	if set(t2t_cols) != set(dict_samples.values()):
		raise ValueError('\n\t# The T2T species names are inconsistent with the dictionary values.' \
						 '\n\t# Perhaps the T2T columns are not selected properly.')
	
	# Concatenate values for T2T cols.
	df['Label'] = df[t2t_cols].astype(str).agg(''.join, axis=1)

	# Label for species-specific NUMTs (ss).
	label_ss = ''.join([ "1" if x == species else "0" for x in t2t_cols ])
	if not quiet:
		print( f"# The species-specific label for '{assembly}' ({species}) NUMTs is:", label_ss)

	# Use label to SS NUMTs.
	df_ss = df[df['Label']==label_ss]#.drop('Label', axis=1)
	if not quiet:
		print( f"# Number of '{assembly}' NUMTs that are species-specific:", df_ss.shape[0])
		print('\n')
	
	# Get the subset of columns for non-T2T population.
	subset_cols = list( set(df_ss.columns) - set(dict_samples.values()) - set(['NUMT_ID','Label']) )

	# Unit test: are these the correct number of non-T2T sample columns?
	if len(subset_cols) != len(df.columns) - len(t2t_cols) - len(['NUMT_ID', 'Label']):
		# Expected number of columns
		print( len(subset_cols), len(t2t_cols), len(['NUMT_ID', 'Label']), len(df.columns) )
		raise ValueError("\n\t# The number of sample columns is different than expected."\
				   		 f"\n\t# Verify the expression: {len(subset_cols)} + {len(t2t_cols)} + {len(['NUMT_ID', 'Label'])} = {len(df.columns)}")
	
	# Concatenate values for non-T2T population cols.
	df_ss = df_ss.copy()
	df_ss['Values'] = df_ss[subset_cols].astype(str).agg(''.join, axis=1)
	df_ss['Classification'] = df_ss['Values'].apply(classify_values)

	# Prepare output.
	df_ss['Assembly'] = assembly
	df_ss['Species'] = 'S. orangutan' if species == 'Sorang' else 'B. orangutan' if species == 'Borang' else species
	df_ss['Total_samples'] = len(subset_cols) + 1 	# Includes T2T
	df_values = df_ss[['Assembly', 'Species', 'NUMT_ID', 'Label', 'Values', 'Total_samples', 'Classification']]

	return df_values
	

# Run for all species.
list_df = []
for sample in dict_samples.keys():
	# Get SS NUMTs.
	print(f"# Getting species-specific NUMTs for {sample}")
	df = get_species_specific(sample, quiet=True)
	
	list_df.append(df)
df = pd.concat(list_df)

# Output values.
file_out = "values.species_specific.tsv"
df.to_csv( file_out, index=None, sep='\t' )
print(f"# Succesfully output to {file_out}")

# Create summary of values.
summary = df.groupby('Species')['Classification'].value_counts().reset_index(name='Count')
summary['Total_samples'] = summary['Species'].map(df.groupby('Species')['Total_samples'].first())

# Specify order to match tree.
species_order = list(df['Species'].unique())
summary['Species'] = pd.Categorical(summary['Species'], categories=species_order, ordered=True)
summary = summary.sort_values(['Species', 'Classification']).reset_index(drop=True)

# Output summary of values.
file_out = "summary.species_specific.tsv"
summary.to_csv( file_out, index=None, sep='\t' )
print(f"# Succesfully output to {file_out}")