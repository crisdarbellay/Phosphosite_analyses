import os
import shutil
import pandas as pd

# Define paths
selection_path = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/to_test/selection_full.txt"
source_dir = "/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/fastas-ball-30"
target_dir = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/to_test/fastas_ball_30"

# Read the data file with a regex separator for whitespace
df = pd.read_csv(selection_path, sep=r'\s+')

# Print the column names to verify
print("Column names in the file:", df.columns)

# Strip whitespace from column names
df.columns = df.columns.str.strip()

# Verify column names again
print("Trimmed column names in the file:", df.columns)

# Define the exact names of the columns based on the structure provided
gene_column = 'Gene'
site_column = 'site'

# Check if the columns exist in the DataFrame
if gene_column not in df.columns or site_column not in df.columns:
    print(f"Available columns: {df.columns}")
    raise KeyError(f"Columns '{gene_column}' and '{site_column}' must be present in the data file")

# Ensure the target directory exists
os.makedirs(target_dir, exist_ok=True)

# Iterate over each row and copy the corresponding FASTA files
for index, row in df.iterrows():
    gene = row[gene_column].strip()
    site = row[site_column]
    fasta_filename = f"{gene}-{site}.fasta"
    source_path = os.path.join(source_dir, fasta_filename)
    target_path = os.path.join(target_dir, fasta_filename)
    
    # Check if the source file exists before copying
    if os.path.exists(source_path):
        shutil.copy(source_path, target_path)
        print(f"Copied: {fasta_filename}")
    else:
        print(f"File not found: {fasta_filename}")
