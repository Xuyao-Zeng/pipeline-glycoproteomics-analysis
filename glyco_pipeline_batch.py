import pandas as pd
import numpy as np
import re
import os

# Hide the future warning about downcasting
pd.set_option('future.no_silent_downcasting', True)

def process_batch_glycoproteomics(input_file, output_file, glycan_file="verified_glycans.xlsx", master_label_file="labeled_proteins.xlsx", new_label_file="protein_list_for_labeling.xlsx", rsd_cutoff=0.30, num_replicates=3):
    print(f"Loading batch data from {input_file}...")
    
    # 1. Read the raw Byonic output
    df_raw = pd.read_excel(input_file, sheet_name='original output')
    
    # 2. Initial QC Filtering
    print("Filtering high confidence and removing blanks...")
    df_qc = df_raw.copy()
    df_qc = df_qc[df_qc['Confidence'] == 'High']
    df_qc = df_qc[df_qc['Quan Info'].isna() | (df_qc['Quan Info'] == '') | (df_qc['Quan Info'] == 'No Quan Values')]
    df_qc = df_qc.dropna(subset=['Glycan Composition'])
    
    # ---------------------------------------------------------
    # --- VERIFIED GLYCAN FILTERING (FIXED) ---
    # ---------------------------------------------------------
    valid_glycans = [] 
    if os.path.exists(glycan_file):
        print(f"Applying verified glycan filter from '{glycan_file}'...")
        df_glycans = pd.read_excel(glycan_file)
        
        # SMART HEADER DETECTION: Find the right column even if the name changes
        if 'Glycan Composition' in df_glycans.columns:
            target_col = 'Glycan Composition'
        elif 'glycan' in df_glycans.columns:
            target_col = 'glycan'
        else:
            target_col = df_glycans.columns[0] # Fallback to the very first column
            
        print(f"Reading verified glycans from column: '{target_col}'...")
        valid_glycans = df_glycans[target_col].dropna().drop_duplicates().tolist()
        
        starting_rows = len(df_qc)
        df_qc = df_qc[df_qc['Glycan Composition'].isin(valid_glycans)]
        ending_rows = len(df_qc)
        print(f"✅ Removed {starting_rows - ending_rows} entries with unverified glycans.")
    else:
        print(f"ℹ️ Verified glycans file '{glycan_file}' not found. Skipping glycan filtering.")
    # ---------------------------------------------------------

    # Extract sample abundance columns using the "Bookend" method
    start_col = 'Theo. MH+ [Da]'
    end_col = 'Quan Info'
    
    start_idx = df_qc.columns.get_loc(start_col) + 1
    end_idx = df_qc.columns.get_loc(end_col)
    abundance_cols = list(df_qc.columns[start_idx:end_idx])
    
    if len(abundance_cols) % num_replicates != 0:
        print(f"⚠️ Warning: Found {len(abundance_cols)} abundance columns, which is not a multiple of {num_replicates}!")
        
    # 3. Batch RSD Filtering
    print("Calculating RSD per sample and masking failed triplicates...")
    df_rsd = df_qc.copy()
    
    for i in range(0, len(abundance_cols), num_replicates):
        triplet = abundance_cols[i:i+num_replicates]
        std = df_rsd[triplet].std(axis=1)
        mean = df_rsd[triplet].mean(axis=1)
        rsd = std / mean
        blanks = df_rsd[triplet].isna().sum(axis=1)
        
        fail_mask = ~((rsd <= rsd_cutoff) & (blanks == 0))
        df_rsd.loc[fail_mask, triplet] = np.nan

    df_rsd = df_rsd.dropna(subset=abundance_cols, how='all')
    
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

    if os.path.exists(new_label_file):
        df_new = pd.read_excel(new_label_file)
        if df_new['Acceptable'].isna().any() or (df_new['Acceptable'] == "").any():
            print("\n" + "="*50)
            print(f"⏸️ PIPELINE PAUSED: '{new_label_file}' still has blank labels!")
            print("Please fill out the remaining 1s and 0s, save, and hit Play again.")
            print("="*50)
            return
        else:
            if os.path.exists(master_label_file):
                master_df = pd.read_excel(master_label_file)
            master_df = pd.concat([master_df, df_new], ignore_index=True).drop_duplicates(subset=['Master Protein Accessions'], keep='last')
            master_df.to_excel(master_label_file, index=False)
            os.remove(new_label_file)
            print(f"✅ Automatically merged your new labels into the master '{master_label_file}'!")

    if os.path.exists(master_label_file):
        master_df = pd.read_excel(master_label_file)

    if not master_df.empty:
        labeled_proteins = master_df['Master Protein Accessions'].tolist()
        missing_proteins = current_proteins[~current_proteins['Master Protein Accessions'].isin(labeled_proteins)].copy()
    else:
        missing_proteins = current_proteins.copy()

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
        
    print(f"\n✅ All proteins are labeled! Applying your GO filters...")
    valid_proteins = master_df[master_df['Acceptable'] == 1]['Master Protein Accessions'].tolist()
    df_gsite = df_gsite[df_gsite['Master Protein Accessions'].isin(valid_proteins)]

    # 6. Normalization
    print("\nNormalizing batch abundance data...")
    df_norm = df_gsite.copy()
    norm_cols = []
    
    for col in abundance_cols:
        col_sum = df_norm[col].sum(skipna=True) 
        norm_col_name = col.replace('Abundance', 'Normalized')
        df_norm[norm_col_name] = (df_norm[col] / col_sum) * 100
        norm_cols.append(norm_col_name)
        
    # Calculate Average Intensity per Sample
    avg_sample_cols = []
    for i in range(0, len(abundance_cols), num_replicates):
        triplet = abundance_cols[i:i+num_replicates]
        norm_triplet = [c.replace('Abundance', 'Normalized') for c in triplet]
        sample_num = (i // num_replicates) + 1
        col_name = f'Average Normalized Sample {sample_num}'
        df_norm[col_name] = df_norm[norm_triplet].mean(axis=1)
        avg_sample_cols.append(col_name)
        
    # ---------------------------------------------------------
    # --- PROTEIN SUMMARY & TOP 10 UNION ---
    # ---------------------------------------------------------
    print("Calculating Glycoprotein summaries and extracting Top 10 lists...")
    df_prot_sum = df_norm.groupby('Master Protein Accessions')[norm_cols].sum(min_count=1).reset_index()
    top_proteins_set = set() 
    
    mean_cols = []
    for i in range(0, len(norm_cols), num_replicates):
        triplet = norm_cols[i:i+num_replicates]
        sample_num = (i // num_replicates) + 1
        mean_col = f'Sample {sample_num} Mean'
        std_col = f'Sample {sample_num} Std'
        mean_cols.append(mean_col)
        
        df_prot_sum[mean_col] = df_prot_sum[triplet].mean(axis=1)
        df_prot_sum[std_col] = df_prot_sum[triplet].std(axis=1)
        top10_for_sample = df_prot_sum.nlargest(10, mean_col)['Master Protein Accessions'].tolist()
        top_proteins_set.update(top10_for_sample)
        
    df_top_union = df_prot_sum[df_prot_sum['Master Protein Accessions'].isin(top_proteins_set)].copy()

    # ---------------------------------------------------------
    # --- ORIGINLAB BUBBLE PLOT MATRIX ---
    # ---------------------------------------------------------
    print("Generating flattened matrix for OriginLab Bubble Plot...")
    
    # 1. Determine Y-axis (Glycans) - This will now accurately retain your custom list order!
    glycans_for_plot = valid_glycans if valid_glycans else sorted(df_norm['Glycan Composition'].unique().tolist())
    
    # 2. Determine X-axis (Proteins)
    df_top_union['Overall Mean'] = df_top_union[mean_cols].mean(axis=1)
    sorted_top_df = df_top_union.sort_values('Overall Mean', ascending=False)
    proteins_for_plot = sorted_top_df['Master Protein Accessions'].tolist()
    df_top_union = df_top_union.drop(columns=['Overall Mean'])
    
    # 3. Create the Cartesian Product
    bubble_data = []
    for x_idx, prot in enumerate(proteins_for_plot, start=1):
        for y_idx, glycan in enumerate(glycans_for_plot, start=1):
            bubble_data.append({
                'X (Protein)': x_idx,
                'Y (Glycan)': y_idx,
                'Protein Accession': prot,
                'Glycan Composition': glycan
            })
    df_bubble = pd.DataFrame(bubble_data)
    
    # 4. Calculate Intensities
    df_pg_sum = df_norm.groupby(['Master Protein Accessions', 'Glycan Composition'])[avg_sample_cols].sum().reset_index()
    
    rename_dict = {col: col.replace('Average Normalized ', '') + ' Intensity' for col in avg_sample_cols}
    df_pg_sum = df_pg_sum.rename(columns=rename_dict)
    
    # 5. Merge the intensities into the flat matrix
    df_bubble = pd.merge(df_bubble, df_pg_sum, how='left', 
                         left_on=['Protein Accession', 'Glycan Composition'], 
                         right_on=['Master Protein Accessions', 'Glycan Composition'])
    
    # Force the sorting order exactly by X and Y!
    df_bubble = df_bubble.sort_values(by=['X (Protein)', 'Y (Glycan)']).reset_index(drop=True)
    
    # Clean up the output
    df_bubble = df_bubble.drop(columns=['Master Protein Accessions']).fillna(0)

    # ---------------------------------------------------------
    # 8. Write to Excel Workbook
    print(f"Saving finalized sheets to {output_file}...")
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        df_raw.to_excel(writer, sheet_name='original output', index=False)
        df_norm.to_excel(writer, sheet_name='normalized_average', index=False)
        df_prot_sum.to_excel(writer, sheet_name='glycoprotein_summary', index=False)
        df_top_union.to_excel(writer, sheet_name='top_glycoproteins', index=False)
        df_bubble.to_excel(writer, sheet_name='origin_bubble_plot', index=False)
        
    print("Batch Pipeline fully complete! 🎉")

if __name__ == "__main__":
    INPUT_FILE = "multiple_data.xlsx" # Change to your new batch file name
    OUTPUT_FILE = "Batch_Pipeline_Output.xlsx"
    
    process_batch_glycoproteomics(INPUT_FILE, OUTPUT_FILE)