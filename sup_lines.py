"""
import os

# Définir le dossier à parcourir

directory= "/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test/results"

def remove_lines_with_none(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Filtrer les lignes ne contenant pas "None"
    filtered_lines = [line for line in lines if 'HEADER' not in line]

    with open(file_path, 'w') as file:
        file.writelines(filtered_lines)
for root, dirs, files in os.walk(directory):
    for filename in files:
        if ".pdb" not in filename:
            continue

        file_path = os.path.join(root, filename)
        remove_lines_with_none(file_path)

print("Suppression des lignes contenant 'None' terminée.")
"""
import os

def process_fasta(input_fasta, site, output_fasta):
    with open(input_fasta, 'r') as f_in:
        lines = f_in.readlines()
        
        # Copier la première ligne (header)
        header = lines[0]
        sequence = lines[1].strip()
        
        # Calculer les indices de coupe
        a = max(site - 30, 0)
        b = min(site + 30, len(sequence) - 1)
        new_sequence = sequence[a:b]
        
        with open(output_fasta, 'w') as f_out:
            f_out.write(header)
            f_out.write(new_sequence + '\n')

input_file = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_filtered.txt"
fasta_folder = "/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/fastas-small/"
output_folder = "/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/fastas-ball-30/"

# Assurez-vous que le dossier de sortie existe
os.makedirs(output_folder, exist_ok=True)

with open(input_file, 'r') as f_in:
    lines = f_in.readlines()
    
    for line in lines[1:]:  # Skip header line
        elements = line.strip().split('\t')
        gene = elements[0]
        site = int(elements[1])
        
        input_fasta_files = [
            f"{gene}.fasta",
            f"{gene}-D{site}.fasta",
            f"{gene}-E{site}.fasta"
        ]
        output_fasta_path = [
            f"{gene}-{site}.fasta",
            f"{gene}-D{site}.fasta",
            f"{gene}-E{site}.fasta"
        ]
        for input_fasta_file, output_fasta_file in zip(input_fasta_files, output_fasta_path):
            input_fasta_path = os.path.join(fasta_folder, input_fasta_file)
            output_fasta_path = os.path.join(output_folder, output_fasta_file)
            
            if os.path.exists(input_fasta_path):
                process_fasta(input_fasta_path, site, output_fasta_path)
            else:
                print(f"File {input_fasta_path} does not exist.")
