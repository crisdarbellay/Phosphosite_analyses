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
    """
    sequence_list = list(sequence)  # Convert sequence to a list of characters
    if sequence_list[site] == site_aa:
        # Wild type sequence is the original
        wild_type_sequence = ''.join(sequence_list)

        # Create D mutation
        sequence_list[site] = 'D'
        modified_sequence_d = ''.join(sequence_list)

        # Create E mutation
        sequence_list[site] = 'E'
        modified_sequence_e = ''.join(sequence_list)
    else:
        # If the site does not match, return empty sequences
        wild_type_sequence = ""
        modified_sequence_d = ""
        modified_sequence_e = ""

    return wild_type_sequence, modified_sequence_d, modified_sequence_e

def create_fastas_and_generate_too_test(filtered_results_file, cif_root_folder, output_folder, too_test_file):
    """
    Read the filtered_results file and create FASTA files. 
    Also, write a copy of each successful entry to the too_test file.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    with open(filtered_results_file, 'r') as file, \
         open(too_test_file, 'w', newline='') as too_test:
        
        reader = csv.reader(file, delimiter='\t')
        writer = csv.writer(too_test, delimiter='\t')

        # Read the header and write it to the too_test file
        header = next(reader)
        writer.writerow(header)
        
        for row_num, row in enumerate(reader, start=2):  # start=2 to account for the header row
            Kinase, gene, site, Score, percentage, alpha, beta, iso_B, alpha3, alphaI, hydturn, AFConf, NextAA, orga, site_aa, id_alphafold, in_vivo, in_vitro, *rest = row

            try:
                cif_folder = os.path.join(cif_root_folder, "human_database")
                cif_file_path = os.path.join(cif_folder, f"AF-{id_alphafold}-F1-model_v4.cif.gz")
                
                if not os.path.exists(cif_file_path):
                    raise FileNotFoundError(f"CIF file not found: {cif_file_path}")

                with gzip.open(cif_file_path, 'rt') as cif_file:
                    protein_sequence = extract_protein_sequence_from_cif(cif_file)
                    
                    if len(protein_sequence) == 0:
                        raise ValueError(f"No protein sequence extracted from CIF file: {cif_file_path}")
                    
                    site_idx = int(site) - 1
                    start = max(site_idx - 29, 0)
                    end = min(site_idx + 30, len(protein_sequence))
                    protein_sequence_segment = protein_sequence[start:end]

                    # Adjust the site index to the new sequence segment
                    adjusted_site_idx = site_idx - start

                    # Perform phosphomimetic simulation
                    wild_type_sequence, modified_sequence_d, modified_sequence_e = phosphomimetic(protein_sequence_segment, adjusted_site_idx, site_aa)

                    if not wild_type_sequence:
                        raise ValueError(f"Site amino acid does not match the expected residue at site {site} in {gene}.")

                    fasta_filename_base = f"{gene}-{site}"
                    fasta_d_name = f"{gene}-D{site}"
                    fasta_e_name = f"{gene}-E{site}"
                    # Write wild type sequence
                    with open(os.path.join(output_folder, f"{fasta_filename_base}.fasta"), 'w') as fasta_file:
                        fasta_file.write(f">{fasta_filename_base}-WT\n{wild_type_sequence}\n")
                    
                    # Write D mutant sequence
                    with open(os.path.join(output_folder, f"{fasta_d_name}.fasta"), 'w') as fasta_file_d:
                        fasta_file_d.write(f">{fasta_d_name}\n{modified_sequence_d}\n")
                    
                    # Write E mutant sequence
                    with open(os.path.join(output_folder, f"{fasta_e_name}.fasta"), 'w') as fasta_file_e:
                        fasta_file_e.write(f">{fasta_e_name}\n{modified_sequence_e}\n")

                    # Write the successful row to the too_test file
                    writer.writerow(row)
                    
            except Exception as e:
                print(f"Error processing row {row_num} ({gene}-{site}): {e}")
                continue

# Example usage
cif_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign"
output_folder = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/3_algorithm/little_fastas"
filtered_results_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/3_algorithm/results/filtered_results.txt"
too_test_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/3_algorithm/results/too_test.txt"

create_fastas_and_generate_too_test(filtered_results_file, cif_folder, output_folder, too_test_file)

print("FASTA files creation and too_test file generation completed.")