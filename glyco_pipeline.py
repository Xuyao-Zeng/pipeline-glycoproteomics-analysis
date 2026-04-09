import pandas as pd
import numpy as np
import re
import os

def process_glycoproteomics_data(input_file, output_file, master_label_file="labeled_proteins.xlsx", new_label_file="protein_list_for_labeling.xlsx", rsd_cutoff=0.30):
    print(f"Loading raw data from {input_file}...")
    
    # 1. Read the raw Byonic output
    df_raw = pd.read_excel(input_file, sheet_name='original output')
    
    # 2. Initial QC Filtering
    print("Filtering high confidence and removing blanks...")
    df_qc = df_raw.copy()
    df_qc = df_qc[df_qc['Confidence'] == 'High']
    df_qc = df_qc[df_qc['Quan Info'].isna() | (df_qc['Quan Info'] == '') | (df_qc['Quan Info'] == 'No Quan Values')]
    df_qc = df_qc.dropna(subset=['Glycan Composition'])
    
    # Extract abundance columns
    abundance_cols = [col for col in df_qc.columns if 'Abundance' in col]
    
    # 3. RSD Filtering
    print("Calculating RSD and filtering...")
    df_rsd = df_qc.copy()
    df_rsd['std'] = df_rsd[abundance_cols].std(axis=1)
    df_rsd['mean'] = df_rsd[abundance_cols].mean(axis=1)
    df_rsd['rsd'] = df_rsd['std'] / df_rsd['mean']
    df_rsd['blank'] = df_rsd[abundance_cols].isna().sum(axis=1)
    df_rsd = df_rsd[(df_rsd['rsd'] <= rsd_cutoff) & (df_rsd['blank'] == 0)]
    
    cols_to_keep = ['Annotated Sequence', 'Master Protein Accessions'] + abundance_cols + ['Glycan Composition', 'Position in Protein']
    df_rsd = df_rsd[cols_to_keep]
    
    # 4. Glycan Data Conversion
    print("Parsing Glycan Compositions...")
    df_gconvert = df_rsd.copy()
    residues = ['Hex', 'HexNAc', 'NeuAc', 'NeuGc', 'Fuc']
    for res in residues:
        pattern = rf'{res}\((\d+)\)'
        df_gconvert[res] = df_gconvert['Glycan Composition'].str.extract(pattern).fillna(0).astype(int)
    
    # 5. Glycosylation Site Mapping
    print("Calculating true glycosylation sites...")
    df_gsite = df_gconvert.copy()
    
    def find_glycosylation_site(row):
        start_pos = row['Position in Protein']
        seq_raw = str(row['Annotated Sequence'])
        if '.' in seq_raw:
            clean_seq = seq_raw.split('.')[1]
        else:
            clean_seq = seq_raw
            
        match = re.search(r'N[^P]([ST]|$)', clean_seq)
        if match:
            return start_pos + match.start()
        return np.nan

    df_gsite['Gsite'] = df_gsite.apply(find_glycosylation_site, axis=1)
    
    # ---------------------------------------------------------
    # --- SMART DATABASE LABELING LOGIC ---
    # ---------------------------------------------------------
    current_proteins = pd.DataFrame(df_gsite['Master Protein Accessions'].unique(), columns=['Master Protein Accessions'])
    master_df = pd.DataFrame(columns=['Master Protein Accessions', 'Acceptable'])

    # Step A: Check if you just finished labeling a new batch
    if os.path.exists(new_label_file):
        df_new = pd.read_excel(new_label_file)
        # Check if there are still blanks
        if df_new['Acceptable'].isna().any() or (df_new['Acceptable'] == "").any():
            print("\n" + "="*50)
            print(f"⏸️ PIPELINE PAUSED: '{new_label_file}' still has blank labels!")
            print("Please fill out the remaining 1s and 0s, save, and hit Play again.")
            print("="*50)
            return
        else:
            # Load master, append the new ones, and save to master
            if os.path.exists(master_label_file):
                master_df = pd.read_excel(master_label_file)
            
            master_df = pd.concat([master_df, df_new], ignore_index=True).drop_duplicates(subset=['Master Protein Accessions'], keep='last')
            master_df.to_excel(master_label_file, index=False)
            os.remove(new_label_file) # Clean up the temporary file!
            print(f"✅ Automatically merged your new labels into the master '{master_label_file}'!")

    # Step B: Load the master list (if it exists)
    if os.path.exists(master_label_file):
        master_df = pd.read_excel(master_label_file)

    # Step C: Find proteins in the CURRENT dataset that are NOT in the master list
    if not master_df.empty:
        labeled_proteins = master_df['Master Protein Accessions'].tolist()
        missing_proteins = current_proteins[~current_proteins['Master Protein Accessions'].isin(labeled_proteins)].copy()
    else:
        missing_proteins = current_proteins.copy() # If no master file exists, all are missing

    # Step D: If there are missing proteins, spit them out and pause!
    if not missing_proteins.empty:
        missing_proteins['Acceptable'] = ""
        missing_proteins.to_excel(new_label_file, index=False)
        print("\n" + "="*50)
        print("⏸️  PIPELINE PAUSED: NEW PROTEINS DETECTED")
        print("="*50)
        print(f"I found {len(missing_proteins)} new proteins that aren't in your master list.")
        print(f"1. I have exported ONLY these new proteins to '{new_label_file}'.")
        print(f"2. Please open it, type '1' or '0' in the 'Acceptable' column, and save it.")
        print(f"3. Hit Play again. I will handle merging it into your master list automatically!")
        return
        
    # Step E: If we reach here, ALL proteins are safely labeled!
    print(f"\n✅ All proteins are labeled! Applying your GO filters...")
    valid_proteins = master_df[master_df['Acceptable'] == 1]['Master Protein Accessions'].tolist()
    
    starting_rows = len(df_gsite)
    df_gsite = df_gsite[df_gsite['Master Protein Accessions'].isin(valid_proteins)]
    ending_rows = len(df_gsite)
    print(f"Removed {starting_rows - ending_rows} misidentified glycopeptides. Proceeding to normalization with {ending_rows} rows.")
    # ---------------------------------------------------------

    # 6. Normalization
    print("\nNormalizing abundance data...")
    df_norm = df_gsite.copy()
    norm_cols = []
    
    for col in abundance_cols:
        col_sum = df_norm[col].sum()
        norm_col_name = col.replace('Abundance', 'Normalized')
        df_norm[norm_col_name] = (df_norm[col] / col_sum) * 100
        norm_cols.append(norm_col_name)
        
    # 7. Calculate Average Intensity
    df_norm['Average Normalized'] = df_norm[norm_cols].mean(axis=1)
    
    # 8. Write to Excel Workbook
    print(f"Saving outputs to {output_file}...")
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df_raw.to_excel(writer, sheet_name='original output', index=False)
        df_rsd.to_excel(writer, sheet_name='rsd filtering', index=False)
        df_gconvert.to_excel(writer, sheet_name='gconvert', index=False)
        df_gsite.to_excel(writer, sheet_name='gsite', index=False)
        
        df_norm.to_excel(writer, sheet_name='normalized_average', index=False)
        
        df_triplicates = df_norm.drop(columns=['Average Normalized'])
        df_triplicates.to_excel(writer, sheet_name='normalized_triplicates', index=False)
        
    print("Pipeline fully complete! 🎉")

if __name__ == "__main__":
    INPUT_FILE = "demonstration10172025.xlsx"
    OUTPUT_FILE = "Pipeline_Debug_Output.xlsx"
    
    process_glycoproteomics_data(INPUT_FILE, OUTPUT_FILE)
