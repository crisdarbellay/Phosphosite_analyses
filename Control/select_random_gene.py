#This code select random gene and site to test them using my algorithm
import os
import gzip
import csv
import random

def extract_protein_sequence_from_cif(cif_file_path):
    """
    Extract the protein sequence from a CIF file.
    """
    protein_sequence = ""
    with gzip.open(cif_file_path, 'rt') as cif_file:
        for line in cif_file:
            if line.startswith("ATOM"):
                parts = line.split()
                if len(parts) >= 6 and parts[3] == 'CA':
                    amino_acid_code = parts[-1]
                    protein_sequence += amino_acid_code
    return protein_sequence

def create_random_site(protein_sequence, exclude_site):
    """
    Create a random site for modification, excluding the given site.
    """
    exclude_index = int(exclude_site[1:]) - 1
    available_sites = [i for i in range(len(protein_sequence)) if i != exclude_index]
    random_site_index = random.choice(available_sites)
    random_site = protein_sequence[random_site_index] + str(random_site_index + 1)
    return random_site

def process_genes(input_file, output_file, cif_root_folder):
    """
    Process genes, select random genes, and generate new sequences with modified sites.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        header = next(reader)
        writer.writerow(header)

        human_genes = [row for row in reader if row[3] == 'human' and row[8] == 'human']
        selected_genes = random.sample(human_genes, 5000)

        for row in selected_genes:
            gene_info = dict(zip(header, row))
            sub_acc_id = gene_info['SUB_ACC_ID']
            sub_mod_rsd = gene_info['SUB_MOD_RSD']

            cif_file_path = os.path.join(cif_root_folder, f"AF-{sub_acc_id}-F1-model_v4.cif.gz")
            try:
                protein_sequence = extract_protein_sequence_from_cif(cif_file_path)
                new_site = create_random_site(protein_sequence, sub_mod_rsd)

                gene_info['SUB_MOD_RSD'] = new_site
                writer.writerow(gene_info.values())

            except Exception as e:
                print(f"Error processing {sub_acc_id}: {e}")
                continue

# Example usage
input_file = r"/mnt/c/Users/crisd/Downloads/Kinase_Substrate_Dataset/Kinase_Substrate_Dataset" 
output_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/control2/1_algorithm/controls.txt"
cif_root_folder =  r"/mnt/c/Users/crisd/Desktop/ProteinDesign/human_database"

process_genes(input_file, output_file, cif_root_folder)
print("Processing completed.")
