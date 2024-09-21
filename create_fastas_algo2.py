import os
import gzip
import csv

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

def phosphomimetic(sequence, site, site_aa):
    """
    Perform phosphomimetic simulation on the given sequence at the specified site.

    Parameters:
    - sequence (str): Protein sequence.
    - site (int): Site of phosphorylation.

    Returns:
    - str: Modified protein sequence.
    """
    sequence_list = list(sequence)  # Convert sequence to a list of characters
    site = int(site) - 1  # Convert to 0-based index
    if site_aa == sequence_list[site]:
        sequence_list[site] = 'D'  # Modify the site of phosphorylation to 'D'
        modified_sequence_d = ''.join(sequence_list)  # Convert the list back to a string

        sequence_list[site] = 'E'  # Modify the site of phosphorylation to 'E'
        modified_sequence_e = ''.join(sequence_list)  # Convert the list back to a string
    elif sequence_list[399]==site_aa:
        sequence_list[399] = 'D'  # Modify the site of phosphorylation to 'D'
        modified_sequence_d = ''.join(sequence_list)  # Convert the list back to a string

        sequence_list[399] = 'E'  # Modify the site of phosphorylation to 'E'
        modified_sequence_e = ''.join(sequence_list)  # Convert the list back to a string   
    else:     
        modified_sequence_d = ""
        modified_sequence_e = ""

    return modified_sequence_d, modified_sequence_e

def create_fastas(missing_predictions_file, cif_root_folder, output_folder, to_test_file, exceptions_file,too_long):
    """
    Process the missing predictions file to generate FASTA files and phosphomimetics.

    Parameters:
    - missing_predictions_file (str): Path to the missing predictions file.
    - cif_folder (str): Path to the folder containing CIF files.
    - output_folder (str): Path to the folder to save the generated FASTA files.
    - to_test_file (str): Path to the file where successfully created FASTA entries will be logged.
    - exceptions_file (str): Path to the file where exceptions will be logged.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    used_genes = set()  # Keep track of used genes
    total_gene = set()

    with open(missing_predictions_file, 'r') as file, \
         open(exceptions_file, 'w', newline='') as exceptions, \
         open(to_test_file, 'w', newline='') as to_test, \
         open(too_long, 'w', newline='') as too_long_file:

        reader = csv.reader(file, delimiter='\t')
        header = next(reader)  # Skip header

        exceptions_writer = csv.writer(exceptions, delimiter='\t')
        to_test_writer = csv.writer(to_test, delimiter='\t')
        too_long_writer = csv.writer(too_long_file, delimiter='\t')

        exceptions_writer.writerow(header)  # Write the same header for exceptions
        to_test_writer.writerow(header)  # Write the same header for to_test
        too_long_writer.writerow(header)  # Write the same header for too_long

        for row in reader:
            Kinase, gene, site, Score,percentage, alpha, beta, iso_B, alpha3, alphaI, hydturn, AFConf, NextAA, orga, site_aa, id_alphafold, in_vivo, in_vitro = row

            try:
                cif_folder = os.path.join(cif_root_folder, "human_database")
                cif_file_path = os.path.join(cif_folder, f"AF-{id_alphafold}-F1-model_v4.cif.gz")
                with gzip.open(cif_file_path, 'rt') as cif_file:
                    protein_sequence = extract_protein_sequence_from_cif(cif_file)
                    if len(protein_sequence) > 800:
                        too_long_writer.writerow(row)
                        site_idx = int(site) - 1
                        start = max(site_idx - 399, 0)
                        end = min(site_idx + 400, len(protein_sequence)-1)
                        protein_sequence = protein_sequence[start:end]
                        fasta_filename = f"{gene}-{site}"
                    else:
                        fasta_filename = gene
                    if gene not in total_gene:
                        total_gene.add(gene)
            except Exception as e:
                exceptions_writer.writerow(row)
                continue

            # Perform phosphomimetic simulation
            try:
                modified_sequence_d, modified_sequence_e = phosphomimetic(protein_sequence, site, site_aa)
            except Exception as e:
                exceptions_writer.writerow(row)
                continue

            # Write phosphomimetic sequences to FASTA files
            with open(os.path.join(output_folder, f"{gene}-D{site}.fasta"), 'w') as fasta_file_d:
                fasta_file_d.write(f">{gene}-D{site}\n{modified_sequence_d}\n")
            with open(os.path.join(output_folder, f"{gene}-E{site}.fasta"), 'w') as fasta_file_e:
                fasta_file_e.write(f">{gene}-E{site}\n{modified_sequence_e}\n")
            with open(os.path.join(output_folder, f"{fasta_filename}.fasta"), 'w') as fasta_file:
                fasta_file.write(f">{fasta_filename}\n{protein_sequence}\n")
            used_genes.add(gene)

            to_test_writer.writerow(row)

    print(len(used_genes))
    print(len(total_gene))

# Example usage
missing_predictions_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/selection_no_doublon.txt"
cif_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign"
output_folder = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/all_fastas_to_test"
to_test_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/to_test.txt"
exceptions_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/exceptions.txt"
too_long = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/too_long.txt"
create_fastas(missing_predictions_file, cif_folder, output_folder, to_test_file, exceptions_file,too_long)

print("FASTA files creation completed.")
