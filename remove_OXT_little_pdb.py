import os

def remove_last_oxt_line(file_path):
    """
    This function removes the last line of the file if it contains 'OXT'.
    
    :param file_path: Path to the PDB file.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Check if the last line contains 'OXT'
    if lines and 'OXT' in lines[-1]:
        lines = lines[:-1]  # Remove the last line

        # Rewrite the file without the last line
        with open(file_path, 'w') as file:
            file.writelines(lines)

def process_all_pdbs_in_directory(root_directory):
    """
    Process all PDB files in the given directory and its subdirectories to remove the last 'OXT' line if present.
    
    :param root_directory: Root directory containing subdirectories with PDB files.
    """
    for subdir, _, files in os.walk(root_directory):
        for file in files:
            if file.endswith(".pdb") and "_relaxed_rank_001_" in file:
                file_path = os.path.join(subdir, file)
                
                # Process the PDB file to remove the last 'OXT' line if present
                remove_last_oxt_line(file_path)
                
                print(f"Processed file: {file_path}")

# Example usage
root_directory = '/mnt/c/Users/crisd/Desktop/mass_test_final/pdb'
process_all_pdbs_in_directory(root_directory)
