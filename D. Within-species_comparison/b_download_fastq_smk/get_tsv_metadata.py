import pandas as pd
import yaml
import json
import glob
import os
import re

# Settings
INPUT_FILES = glob.glob("src/*.tsv")
OUTPUT_YAML = "config/samples.yaml"

def clean_sra_string(s):
    """Parses JSON lists or raw strings into a clean list of IDs."""
    if pd.isna(s) or s == "":
        return []
    s = str(s).strip()
    # Handle the ["SRR123","SRR456"] or "[""ERR123""]" formats
    if s.startswith("[") and s.endswith("]"):
        try:
            # Normalize double-double quotes
            s_json = s.replace('""', '"')
            return json.loads(s_json)
        except:
            # Fallback: regex find anything looking like an accession
            return re.findall(r'[S|E|D]R[R|X]\d+', s)
    # Handle comma-separated strings (like in s233_extra_download.tsv)
    if "," in s:
        return [x.strip() for x in s.split(",")]
    # Handle single IDs
    return [s]

def combine_and_convert():
    combined_data = {}
    all_srr_ids = set()
    processed_files = []

    for tsv in INPUT_FILES:
        try:
            df = pd.read_csv(tsv, sep='\t')
            processed_files.append(os.path.basename(tsv))
            
            # Normalize column names for this specific file
            cols = {c.lower(): c for c in df.columns}
            
            # Identify ID column (handles 'entity:sample_id' and 'entity:sample')
            id_col = cols.get('entity:sample_id') or cols.get('entity:sample')
            species_col = cols.get('common_name')
            sra_list_col = cols.get('sra_list')
            srx_col = cols.get('srx')
            sex_col = cols.get('sex')

            for _, row in df.iterrows():
                if pd.isna(row[id_col]): continue
                
                indiv_name = str(row[id_col]).strip()
                species = str(row.get(species_col, "Unknown")).strip().replace(" ", "_")
                sex = str(row.get(sex_col, "unknown")).strip().lower()
                
                if indiv_name not in combined_data:
                    combined_data[indiv_name] = {
                        'species': species,
                        'sex': sex if sex != 'nan' else 'unknown',
                        'runs': set()
                    }

                # Extract from SRA_List
                if sra_list_col:
                    for run in clean_sra_string(row[sra_list_col]):
                        combined_data[indiv_name]['runs'].add(run)
                        if "RR" in run: all_srr_ids.add(run)

                # Extract from SRX/ERX column
                if srx_col:
                    for run in clean_sra_string(row[srx_col]):
                        combined_data[indiv_name]['runs'].add(run)
                        # Note: If it's an SRX/ERX, it might represent multiple SRR_1 files
                        if "RR" in run: all_srr_ids.add(run)

        except Exception as e:
            print(f"# Error reading {tsv}: {e}")

    # Prepare final dictionary
    final_output = {'individuals': {}}
    for indiv, data in combined_data.items():
        final_output['individuals'][indiv] = {
            'species': data['species'],
            'sex': data['sex'],
            'runs': sorted(list(data['runs']))
        }

    # Write YAML
    os.makedirs(os.path.dirname(OUTPUT_YAML), exist_ok=True)
    with open(OUTPUT_YAML, 'w') as f:
        yaml.dump(final_output, f, sort_keys=False, default_flow_style=False)

    # FINAL LOGS
    print(f"# Metadata processing complete.")
    print(f"# Files processed: {', '.join(processed_files)}")
    print(f"# Total individuals found: {len(final_output['individuals'])}")
    print(f"# Total unique Run/Experiment IDs found: {sum(len(v['runs']) for v in final_output['individuals'].values())}")

if __name__ == "__main__":
    combine_and_convert()
