import os
import gzip

def extract_protein_sequence_from_cif(cif_file):
    """
    Extract the protein sequence from a CIF file.
    """
    protein_sequence = ""
    amino_acid_code = None

    for line in cif_file:
        if line.startswith("ATOM"):
            parts = line.split()
            if len(parts) >= 6 and parts[3] == 'CA':
                amino_acid_code = parts[-1]
                protein_sequence += amino_acid_code

    return protein_sequence

def phosphomimetic(sequence, site,site_aa):
    """
    Perform phosphomimetic simulation on the given sequence at the specified site.

    Parameters:
    - sequence (str): Protein sequence.
    - site (int): Site of phosphorylation.

    Returns:
    - str: Modified protein sequence.
    """
    sequence_list = list(sequence)  # Convert sequence to a list of characters
    site = int(site)
    if site_aa == sequence_list[site - 1]:
        sequence_list[site - 1] = 'D'  # Modify the site of phosphorylation to 'D'
        modified_sequence_d = ''.join(sequence_list)  # Convert the list back to a string

        sequence_list[site - 1] = 'E'  # Modify the site of phosphorylation to 'E'
        modified_sequence_e = ''.join(sequence_list)  # Convert the list back to a string
    elif site_aa == sequence_list[site]:
        sequence_list[site] = 'D'  # Modify the site of phosphorylation to 'D'
        modified_sequence_d = ''.join(sequence_list)  # Convert the list back to a string

        sequence_list[site] = 'E'  # Modify the site of phosphorylation to 'E'
        modified_sequence_e = ''.join(sequence_list)  # Convert the list back to a string
    elif site_aa == sequence_list[site+1]:
        sequence_list[site+1] = 'D'  # Modify the site of phosphorylation to 'D'
        modified_sequence_d = ''.join(sequence_list)  # Convert the list back to a string

        sequence_list[site+1] = 'E'  # Modify the site of phosphorylation to 'E'
        modified_sequence_e = ''.join(sequence_list)  # Convert the list back to a string
    else:
        modified_sequence_d = ""
        modified_sequence_e = ""


    return modified_sequence_d, modified_sequence_e

def process_phosphosite_file(input_file, cif_folder):
    """
    Process the phosphorylation site file to generate FASTA files and phosphomimetics.

    Parameters:
    - input_file (str): Path to the input phosphorylation site file.
    - cif_folder (str): Path to the folder containing CIF files.
    """
    used_genes = set()  # Keep track of used genes
    total_gene = set()

    with open(input_file, 'r') as file:
        for line in file:
            try:
                gene,site,Score,AFConf,NextAA,in_vivo,in_vitro,id_alphafold,site_aa = line.split()
            except:
                continue
  # Add gene to used genes
            try:
                cif_file_path = os.path.join(cif_folder, f"AF-{id_alphafold}-F1-model_v4.cif.gz")
                with gzip.open(cif_file_path, 'rt') as cif_file:
                    protein_sequence = extract_protein_sequence_from_cif(cif_file)
                    if gene not in total_gene:
                        total_gene.add(gene)
                    if len(protein_sequence)>1000:
                        continue
            except:
                continue

            # Perform phosphomimetic simulation
            try:
                modified_sequence_d, modified_sequence_e = phosphomimetic(protein_sequence, site,site_aa)
            except:
                continue

            # Write phosphomimetic sequences to FASTA files
            with open(f"/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/fastas-small/{gene}-D{site}.fasta", 'w') as fasta_file_d:
                fasta_file_d.write(f">{gene}-D{site}\n{modified_sequence_d}\n")
            with open(f"/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/fastas-small/{gene}-E{site}.fasta", 'w') as fasta_file_e:
                fasta_file_e.write(f">{gene}-E{site}\n{modified_sequence_e}\n")

            if gene not in used_genes:  # Skip if gene has already been used
                with open(f"/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/fastas-small/{gene}.fasta", 'w') as fasta_file:
                    fasta_file.write(f">{gene}\n{protein_sequence}\n")
            used_genes.add(gene)
        print(len(used_genes))
        print(len(total_gene))
# Exemple d'utilisation :
input_file = '/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results.txt'
cif_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/human_database"
process_phosphosite_file(input_file, cif_folder)
