import pandas as pd
import ast
import pysam
from Bio.Seq import Seq
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import sys

def process_row(args):
    row, GenomeDict2, FLANK_SIZE = args
    mykmerSize = 14
    dissimilarities = []
    similarities = []

    contig = str(row['Loci'].split(":")[2])
    start = str(int(row['Loci'].split(":")[3].split("-")[0]) - FLANK_SIZE)
    end = str(int(row['Loci'].split(":")[3].split("-")[1]) + FLANK_SIZE)

    focus_species = str(row['Query_Species'])

    try:
        sequence = ''.join(pysam.faidx(GenomeDict2[focus_species], f"{contig}:{start}-{end}").split()[1:])
    except Exception as e:
        print(f"Error fetching reference sequence: {e}")
        return [], []

    for match in ast.literal_eval(str(row['Matches'])):
        querySpecies = str(row['DB_Species'])

        coordinate = str(match.split("::")[1])

        try:
            if str(match.split("::")[0]) == '+':
                sequence2 = ''.join(pysam.faidx(GenomeDict2[querySpecies], coordinate).split()[1:])
            else:
                presequence = Seq(''.join(pysam.faidx(GenomeDict2[querySpecies], coordinate).split()[1:]))
                sequence2 = str(presequence.reverse_complement())
        except Exception as e:
            print(f"Error fetching query sequence: {e}")
            continue

        kmers = []
        insertionkmers = []
        i = 0
        while i < (len(sequence) - (mykmerSize - 1)):
            kmers.append(str(sequence[i:i + mykmerSize]).upper())
            insertionkmers.append(str(sequence2[i:i + mykmerSize]).upper())
            i += 1

        uniqueKmers = set(kmers + insertionkmers)
        tempDF = pd.DataFrame(0, index=['Reference', 'MEI'], columns=[x for x in uniqueKmers])

        for kmer in kmers:
            tempDF.at['Reference', kmer] += 1
        for kmer2 in insertionkmers:
            tempDF.at['MEI', kmer2] += 1

        tempDF.loc['Diff'] = tempDF.loc['Reference'] - tempDF.loc['MEI']
        tempDF.loc['Sum'] = tempDF.loc['Reference'] + tempDF.loc['MEI']
        dissimilarities.append(sum(abs(tempDF.loc['Diff'])) / sum(abs(tempDF.loc['Sum'])))
        similarities.append(1 - sum(abs(tempDF.loc['Diff'])) / sum(abs(tempDF.loc['Sum'])))

    return dissimilarities, similarities

def expand_columns(file_path, output_file_path):
    df = pd.read_table(file_path)
    df['Query_CommonName'] = [{
        'CHM13':'Human', 'mGorGor1':'Gorilla', 'mPanPan1':'Bonobo', 'mPanTro3':'Chimpanzee', 'mSymSyn1':'Siamang', 'mPonAbe1':'Sorang', 'mPonPyg2':'Borang'
    }[x] for x in df['Query_Species']]
    df['DB_CommonName'] = [{
        'CHM13':'Human', 'mGorGor1':'Gorilla', 'mPanPan1':'Bonobo', 'mPanTro3':'Chimpanzee', 'mSymSyn1':'Siamang', 'mPonAbe1':'Sorang', 'mPonPyg2':'Borang'
    }[x] for x in df['DB_Species']]
    
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
    
    df_expanded.to_csv(output_file_path, sep='\t', index=False)
    return df_expanded

if __name__ == '__main__':
    # Ensure the user has provided a FLANK_SIZE argument
    if len(sys.argv) != 2:
        print("Usage: python script.py <FLANK_SIZE>")
        sys.exit(1)

    # Get the FLANK_SIZE from the command-line arguments
    FLANK_SIZE = int(sys.argv[1])
    DIR=f"results_{FLANK_SIZE}bp/"

    dataDF = pd.read_table(f"{DIR}blat_matches.flanks{FLANK_SIZE}.tab")

    GenomeDict2 = {
        "mGorGor1": "../../data/refs/mGorGor1.pri.cur.20231122.fasta",
        "mPanPan1": "../../data/refs/mPanPan1.pri.cur.20231122.fasta",
        "mPanTro3": "../../data/refs/mPanTro3.pri.cur.20231122.fasta",
        "mPonAbe1": "../../data/refs/mPonAbe1.pri.cur.20231205.fasta",
        "mPonPyg2": "../../data/refs/mPonPyg2.pri.cur.20231122.fasta",
        "mSymSyn1": "../../data/refs/mSymSyn1.pri.cur.20231205.fasta",
        "CHM13": "../../data/refs/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna"
    }

    rows = [(dataDF.loc[i], GenomeDict2, FLANK_SIZE) for i in dataDF.index]

    with Pool(cpu_count()) as pool:
        results = list(tqdm(pool.imap(process_row, rows), total=len(rows), desc="Processing rows", unit="row"))

    seqDiss, seqSims = zip(*results)

    dataDF['Hit_Dissimilarity'] = seqDiss
    dataDF['Hit_Similarity'] = seqSims
    dataDF.to_csv(f"{DIR}matches_similarity_results.flanks{FLANK_SIZE}.tab", index=None, sep='\t')

    expand_columns(file_path=f'{DIR}matches_similarity_results.flanks{FLANK_SIZE}.tab',
                   output_file_path=f'{DIR}expanded_matches_similarity_results.flanks{FLANK_SIZE}.tab')
