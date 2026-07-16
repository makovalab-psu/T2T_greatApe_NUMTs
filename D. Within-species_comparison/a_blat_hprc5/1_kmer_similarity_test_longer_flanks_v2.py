import pandas as pd
import ast
import pysam
import os
from Bio.Seq import Seq
from tqdm import tqdm
import sys

FLANK_SIZE = 500
MY_KMER_SIZE = 14

GenomeDict2 = {
    "mGorGor1": "mGorGor1.pri.cur.20231122.fasta",
    "mPanPan1": "mPanPan1.pri.cur.20231122.fasta",
    "mPanTro3": "mPanTro3.pri.cur.20231122.fasta",
    "mPonAbe1": "mPonAbe1.pri.cur.20231205.fasta",
    "mPonPyg2": "mPonPyg2.pri.cur.20231122.fasta",
    "mSymSyn1": "mSymSyn1.pri.cur.20231205.fasta",
    "CHM13": "GCF_009914755.1_T2T-CHM13v2.0_genomic.fna",
    "HG002": "HG002v1.1.pat.PanSN.fa",
    "HG00597": "HG00597_pat_hprc_r2_v1.0.1.fa",
    "HG01358": "HG01358_pat_hprc_r2_v1.0.1.fa",
    "HG02572": "HG02572_pat_hprc_r2_v1.0.1.fa",
    "HG04184": "HG04184_pat_hprc_r2_v1.0.1.fa",
}

SPECIES_MAP = {
    'CHM13': 'Human', 'mGorGor1': 'Gorilla', 'mPanPan1': 'Bonobo',
    'mPanTro3': 'Chimpanzee', 'mSymSyn1': 'Siamang', 'mPonAbe1': 'Sorang',
    'mPonPyg2': 'Borang'
}

seqDiss = []
seqSims = []

# Read the tab-delimited file into a DataFrame
input_file = f"blat_matches.flanks{FLANK_SIZE}.tab"
dataDF = pd.read_table(input_file)

print(f"Processing {len(dataDF)} loci for k-mer similarity analysis...")

for row in tqdm(dataDF.index, desc="Processing rows", unit="row"):
    dissimilarities = []
    similarities = []

    # Parse out coordinates from 'Loci' column
    loci_str = str(dataDF.at[row, 'Loci'])
    parts = loci_str.split(":")
    
    # Handle complex PanSN naming tokens securely
    contig = ":".join(parts[2:-1]) if len(parts) > 4 else parts[2]
    coord_split = parts[-1].split("-")
    
    # Enforce a 0-based index floor for flanking boundaries
    start = max(0, int(coord_split[0]) - FLANK_SIZE)
    end = int(coord_split[1]) + FLANK_SIZE
    focus_species = str(dataDF.at[row, 'Query_Species'])

    # 1. Fetch Query Reference Sequence cleanly via context manager
    if focus_species not in GenomeDict2:
        seqDiss.append([])
        seqSims.append([])
        continue
        
    try:
        with pysam.FastaFile(GenomeDict2[focus_species]) as fa:
            if contig in fa.references:
                # Safe-crop boundaries against actual contig lengths
                end = min(end, fa.get_reference_length(contig))
                sequence = fa.fetch(contig, start, end).upper()
            else:
                seqDiss.append([])
                seqSims.append([])
                continue
    except Exception as e:
        print(f"Error fetching reference sequence for {loci_str}: {e}")
        raise ValueError("# Skipping!")
        seqDiss.append([])
        seqSims.append([])
        continue

    # 2. Process Target Matches
    db_species = str(dataDF.at[row, 'DB_Species'])
    if db_species not in GenomeDict2:
        raise ValueError("# Skipping!")
        seqDiss.append([])
        seqSims.append([])
        continue

    try:
        matches_list = ast.literal_eval(str(dataDF.at[row, 'Matches']))
    except Exception:
        raise ValueError("# Skipping!")
        seqDiss.append([])
        seqSims.append([])
        continue

    # Open target database file index once per row to maximize memory efficiency
    with pysam.FastaFile(GenomeDict2[db_species]) as fa_db:
        for match in matches_list:
            try:
                match_parts = match.split("::")
                strand = match_parts[0]
                coordinate = match_parts[1]
                
                coord_chrom = ":".join(coordinate.split(":")[:-1]) if coordinate.count(":") > 1 else coordinate.split(":")[0]
                db_start, db_end = map(int, coordinate.split(":")[-1].split("-"))
                
                if coord_chrom not in fa_db.references:
                    continue
                    
                db_start = max(0, db_start)
                db_end = min(db_end, fa_db.get_reference_length(coord_chrom))
                
                seq_raw = fa_db.fetch(coord_chrom, db_start, db_end).upper()
                
                if strand == '+':
                    sequence2 = seq_raw
                else:
                    sequence2 = str(Seq(seq_raw).reverse_complement()).upper()
            except Exception as e:
                print(f"Error fetching database sequence {coordinate}: {e}")
                continue

            # 3. Fast list-comprehension k-mer extraction
            kmers = [sequence[i:i + MY_KMER_SIZE] for i in range(len(sequence) - (MY_KMER_SIZE - 1))]
            insertionkmers = [sequence2[i:i + MY_KMER_SIZE] for i in range(len(sequence2) - (MY_KMER_SIZE - 1))]
            
            unique_kmers = set(kmers + insertionkmers)
            
            if not unique_kmers:
                dissimilarities.append(0.0)
                similarities.append(1.0)
                continue
                
            # Compute similarity counts using lightweight dictionaries
            ref_counts = {k: 0 for k in unique_kmers}
            mei_counts = {k: 0 for k in unique_kmers}
            
            for k in kmers: ref_counts[k] += 1
            for k in insertionkmers: mei_counts[k] += 1
            
            diff_sum = 0
            total_sum = 0
            for k in unique_kmers:
                diff_sum += abs(ref_counts[k] - mei_counts[k])
                total_sum += abs(ref_counts[k] + mei_counts[k])
                
            if total_sum == 0:
                dissimilarities.append(0.0)
                similarities.append(1.0)
            else:
                diss = diff_sum / total_sum
                dissimilarities.append(diss)
                similarities.append(1.0 - diss)

    seqDiss.append(dissimilarities)
    seqSims.append(similarities)

# Update intermediate outputs
dataDF['Hit_Dissimilarity'] = seqDiss
dataDF['Hit_Similarity'] = seqSims

inter_file = f"matches_similarity_results.flanks{FLANK_SIZE}.tab"
dataDF.to_csv(inter_file, index=None, sep='\t')


def expand_columns(file_path, output_file_path):
    """
    Expands rows of the DataFrame based on the lists in 'Matches', Hit_Dissimilarity and 'Hit_Similarity' columns.
    """
    df = pd.read_table(file_path)

    # Clean map fallback for missing common name keys
    df['Query_CommonName'] = df['Query_Species'].map(SPECIES_MAP).fillna(df['Query_Species'])
    df['DB_CommonName'] = df['DB_Species']
    print(df.columns)
    
    df['Matches'] = df['Matches'].apply(ast.literal_eval)
    df['Hit_Dissimilarity'] = df['Hit_Dissimilarity'].apply(ast.literal_eval)
    df['Hit_Similarity'] = df['Hit_Similarity'].apply(ast.literal_eval)

    rows = []
    for i, row in df.iterrows():
        matches = row['Matches']
        dissim = row['Hit_Dissimilarity']
        sim = row['Hit_Similarity']
        
        # Guard rails for anomalous runs or un-calculated rows
        if len(dissim) == 0 or len(sim) == 0:
            continue
            
        if len(matches) != len(dissim):
            print(f"Skipping index {i} due to structural calculation anomalies ({len(matches)} vs {len(dissim)})")
            continue
            
        base_dict = row.to_dict()  # dictionary conversions are 100x faster than row.copy()
        for j in range(len(matches)):
            new_row = base_dict.copy()
            new_row['Matches'] = matches[j]
            new_row['Hit_Dissimilarity'] = dissim[j]
            new_row['Hit_Similarity'] = sim[j]
            new_row['List_Index'] = j
            rows.append(new_row)
            
    df_expanded = pd.DataFrame(rows)
    df_expanded.to_csv(output_file_path, sep='\t', index=False)
    print(f"Expanded mapping file generated successfully: {output_file_path}")
    return df_expanded


# Final Execution Pass
expand_columns(file_path=inter_file,
               output_file_path=f'expanded_matches_similarity_results.flanks{FLANK_SIZE}.tab')

