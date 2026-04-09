import pandas as pd
import numpy as np

def calculate_multi_dim_corr(input_file, output_file):
    print(f"Loading data from {input_file}...")
    
    # Read the data (header=None because the input is a pure numeric matrix)
    df = pd.read_excel(input_file, header=None)
    data_matrix = df.values
    
    ht, wd = data_matrix.shape
    dm = 5 # Number of dimensions (Hex, HexNAc, NeuAc, NeuGc, Fuc)
    
    # Separate the 5 dimension columns from the actual data columns
    dims = data_matrix[:, :dm]
    data_cols = data_matrix[:, dm:]
    
    num_data_cols = wd - dm
    
    # Initialize the output matrix with zeros
    outputtb = np.zeros((num_data_cols, num_data_cols))
    
    print(f"Calculating multi-dimensional correlation for {num_data_cols}x{num_data_cols} matrix...")
    
    for i in range(num_data_cols):
        # Extract the current column and keep it as a 2D column vector for broadcasting
        C1 = data_cols[:, i:i+1]
        
        # Precalculate input1 (element-wise multiplication with all 5 dimensions)
        input1 = dims * C1
        
        # Calculate mean of each dimension column and mean-center the residuals
        avg1 = np.mean(input1, axis=0)
        res1 = input1 - avg1
        
        for j in range(num_data_cols):
            C2 = data_cols[:, j:j+1]
            input2 = dims * C2
            
            avg2 = np.mean(input2, axis=0)
            res2 = input2 - avg2
            
            # Replicating MATLAB's sum(dot(res1, res2))
            # In Python/Numpy, this is the sum of the element-wise products of the matrices
            num = np.sum(res1 * res2)
            den = np.sqrt(np.sum(res1 * res1)) * np.sqrt(np.sum(res2 * res2))
            
            if den == 0:
                outputtb[i, j] = np.nan # Avoid division by zero if vectors are perfectly flat
            else:
                outputtb[i, j] = num / den
                
    # Save the output table
    print(f"Saving output to {output_file}...")
    pd.DataFrame(outputtb).to_excel(output_file, index=False, header=False)
    print("Complete!")

if __name__ == "__main__":
    INPUT_FILE = "Multi_di.xlsx"
    OUTPUT_FILE = "Result.xlsx"
    
    calculate_multi_dim_corr(INPUT_FILE, OUTPUT_FILE)