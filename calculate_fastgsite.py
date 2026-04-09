import pandas as pd

def calculate_gsite(input_file, output_file):
    print(f"Loading data from {input_file}...")
    df = pd.read_excel(input_file)
    
    # Define a function to calculate the exact position
    def find_j_position(row):
        peptide = str(row['Peptide'])
        pro_site = int(row['ProSite'])
        
        # .find() gets the 0-based index of 'J' in the string
        # e.g., if J is the first letter, index is 0. 
        # Adding this to the starting position gives the exact protein site!
        j_index = peptide.find('J')
        
        if j_index != -1:
            return pro_site + j_index
        else:
            return None # Just in case a row is missing the 'J'

    print("Calculating absolute glycosylation sites...")
    df['Gsite'] = df.apply(find_j_position, axis=1)
    
    print(f"Saving output to {output_file}...")
    df.to_excel(output_file, index=False)
    print("Complete! 🎉")

if __name__ == "__main__":
    INPUT_FILE = "fastgsite.xlsx"
    OUTPUT_FILE = "fastgsite_Output.xlsx"
    
    calculate_gsite(INPUT_FILE, OUTPUT_FILE)