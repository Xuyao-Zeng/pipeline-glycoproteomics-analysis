import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def process_site_correlation(input_file, output_file, min_glycans=0):
    print(f"Loading data from {input_file}...")
    df = pd.read_excel(input_file)
    
    # 1. Dynamically identify sample columns (everything after 'Gsite')
    start_idx = df.columns.get_loc('Gsite') + 1
    sample_cols = df.columns[start_idx:]
    
    print(f"Detected {len(sample_cols)} sample columns.")
    
    # 2. Create a unique Site ID (e.g., 'P01833_499')
    df['Site_ID'] = df['Master Protein Accessions'].astype(str) + '_' + df['Gsite'].astype(str)
    
    # 3. Calculate how many glycans exist per site per sample
    counts_list = []
    for col in sample_cols:
        # A glycan is considered present if its abundance is > 0
        valid_mask = df[col].notna() & (df[col] > 0)
        site_counts = df[valid_mask].groupby('Site_ID').size()
        counts_list.append((col, site_counts))
        
    # ---------------------------------------------------------
    # --- PHASE 1: GENERATE HISTOGRAMS ---
    # ---------------------------------------------------------
    print("Generating histograms of glycan counts per site...")
    num_samples = len(sample_cols)
    cols_plot = min(4, num_samples)
    rows_plot = int(np.ceil(num_samples / cols_plot))
    
    fig, axes = plt.subplots(rows_plot, cols_plot, figsize=(cols_plot*4, rows_plot*3.5))
    if num_samples > 1:
        axes = axes.flatten()
    else:
        axes = [axes]
        
    for i, (col_name, counts) in enumerate(counts_list):
        if len(counts) > 0:
            max_val = int(counts.max())
            axes[i].hist(counts, bins=range(1, max_val + 2), edgecolor='black', alpha=0.7, align='left')
        axes[i].set_title(col_name, fontsize=10)
        axes[i].set_xlabel("Number of Glycans", fontsize=9)
        axes[i].set_ylabel("Number of Sites", fontsize=9)
        axes[i].tick_params(axis='both', which='major', labelsize=8)

    # Hide any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    hist_filename = 'glycan_counts_histograms.png'
    plt.savefig(hist_filename, dpi=300)
    print(f"✅ Saved histograms to '{hist_filename}'.")
    
    # Stop here if the user hasn't set a cutoff yet
    if min_glycans <= 0:
        print("\n" + "="*50)
        print("⏸️ PIPELINE PAUSED: REVIEW HISTOGRAMS")
        print("="*50)
        print(f"I have generated '{hist_filename}'.")
        print("Please look at the distribution (you can check the 'Plots' pane in Spyder) and decide on a reasonable cutoff.")
        print("Then, change 'MIN_GLYCANS = 0' to your chosen number at the bottom of this script and hit Play again!")
        return

    # ---------------------------------------------------------
    # --- PHASE 2: MULTI-DIMENSIONAL CORRELATION MATRIX ---
    # ---------------------------------------------------------
    print(f"\nApplying cutoff: Sites must have >= {min_glycans} glycans.")
    
    # Dictionary to keep track of which sites exist in which sample
    site_existence_tracker = {}
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for col, site_counts in counts_list:
            sheet_name = str(col)[:31] # Excel limits sheet names to 31 chars
            print(f"Processing 5D correlation for sample: {col}...")
            
            # Filter sites that meet the cutoff
            valid_sites = site_counts[site_counts >= min_glycans].index.tolist()
            
            if len(valid_sites) < 2:
                print(f"  ⚠️ Not enough valid sites ({len(valid_sites)}) in {col} to correlate. Skipping.")
                continue
                
            print(f"  > Found {len(valid_sites)} valid sites.")
            
            # Record these valid sites for our final summary sheet
            site_existence_tracker[sheet_name] = valid_sites
            
            # Subset the dataframe for this sample's valid sites
            df_sample = df[df['Site_ID'].isin(valid_sites)].copy()
            
            # Pivot the matrix so columns = Sites, rows = Glycan Compositions
            pivot_df = pd.pivot_table(df_sample, 
                                      values=col, 
                                      index=['Hex', 'HexNAc', 'NeuAc', 'NeuGc', 'Fuc', 'Glycan Composition'], 
                                      columns='Site_ID', 
                                      aggfunc='sum', 
                                      fill_value=0).reset_index()
                                      
            dims = pivot_df[['Hex', 'HexNAc', 'NeuAc', 'NeuGc', 'Fuc']].values
            
            site_columns = [c for c in pivot_df.columns if c not in ['Hex', 'HexNAc', 'NeuAc', 'NeuGc', 'Fuc', 'Glycan Composition']]
            data_cols = pivot_df[site_columns].values
            
            num_sites = len(site_columns)
            outputtb = np.zeros((num_sites, num_sites))
            
            # 5D Dot-Product Correlation Math
            for i in range(num_sites):
                C1 = data_cols[:, i:i+1]
                input1 = dims * C1
                avg1 = np.mean(input1, axis=0)
                res1 = input1 - avg1
                
                for j in range(num_sites):
                    C2 = data_cols[:, j:j+1]
                    input2 = dims * C2
                    avg2 = np.mean(input2, axis=0)
                    res2 = input2 - avg2
                    
                    den = np.sqrt(np.sum(res1 * res1)) * np.sqrt(np.sum(res2 * res2))
                    if den == 0:
                        outputtb[i, j] = np.nan
                    else:
                        outputtb[i, j] = np.sum(res1 * res2) / den
                        
            corr_df = pd.DataFrame(outputtb, index=site_columns, columns=site_columns)
            corr_df.to_excel(writer, sheet_name=sheet_name)
            
        # ---------------------------------------------------------
        # --- PHASE 3: SITE EXISTENCE SUMMARY MATRIX ---
        # ---------------------------------------------------------
        if site_existence_tracker:
            print("\nGenerating Site Existence summary matrix...")
            
            # Extract a unique master list of all sites that passed the cutoff in ANY sample
            all_unique_sites = set()
            for sites in site_existence_tracker.values():
                all_unique_sites.update(sites)
                
            # Build the existence dataframe
            existence_data = {'Site_ID': sorted(list(all_unique_sites))}
            
            # Loop through each sample and assign a 1 if the site exists, or a 0 if it doesn't
            for sample_name, valid_sites_list in site_existence_tracker.items():
                existence_data[sample_name] = [1 if site in valid_sites_list else 0 for site in existence_data['Site_ID']]
                
            df_existence = pd.DataFrame(existence_data)
            
            # Optional: Add a "Total Occurrences" column at the very end to sum them up
            df_existence['Total Occurrences'] = df_existence.drop(columns=['Site_ID']).sum(axis=1)
            
            # Save it to the workbook!
            df_existence.to_excel(writer, sheet_name='Site_Existence', index=False)
            print("✅ 'Site_Existence' sheet added successfully.")
            
    print(f"\n✅ All results fully saved to '{output_file}'! 🎉")


if __name__ == "__main__":
    INPUT_FILE = "site_correl.xlsx"
    OUTPUT_FILE = "Site_Correlation_Results.xlsx"
    
    # >>> CHANGE THIS NUMBER ONCE YOU REVIEW THE HISTOGRAMS <<<
    MIN_GLYCANS = 5  
    
    process_site_correlation(INPUT_FILE, OUTPUT_FILE, min_glycans=MIN_GLYCANS)