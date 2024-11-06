import os
import gzip
import csv

def extract_protein_sequence_from_cif(cif_file):
    """
    Extract the protein sequence from a CIF file.
    """
    protein_sequence = ""
    for line in cif_file:
        if line.startswith("ATOM"):
            parts = line.split()
            if len(parts) >= 6 and parts[3] == 'CA':  # Assuming alpha-carbon (CA) line
                amino_acid_code = parts[-1]  # Amino acid code
                protein_sequence += amino_acid_code
    return protein_sequence

def phosphomimetic(sequence, sites, site_aa):
    """
    Perform phosphomimetic simulation on the given sequence at multiple specified sites.
    
    Parameters:
    - sequence (str): Protein sequence.
    - sites (list): List of sites of phosphorylation.

    Returns:
    - str: Modified protein sequence with E in place of phosphorylation sites.
    """
    sequence_list = list(sequence)  # Convert sequence to a list of characters
    for site in sites:
        site_idx = int(site) - 1  # Convert to 0-based index
        if sequence_list[site_idx] == site_aa:
            sequence_list[site_idx] = 'E'  # Modify the site of phosphorylation to 'E'
    return ''.join(sequence_list)  # Convert the list back to a string

def create_fastas(input_file, cif_root_folder, output_folder, to_test_file, exceptions_file, too_long_file, fasta_log_file):
    """
    Process the input file to generate FASTA files and phosphomimetics.
    
    Parameters:
    - input_file (str): Path to the input file.
    - cif_root_folder (str): Path to the folder containing CIF files.
    - output_folder (str): Path to the folder to save the generated FASTA files.
    - to_test_file (str): Path to the file where successfully created FASTA entries will be logged.
    - exceptions_file (str): Path to the file where exceptions will be logged.
    - too_long_file (str): Path to the file for proteins exceeding length limits.
    - fasta_log_file (str): Path to the file that logs the genes and sites of created FASTAs.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    used_genes = set()  # Keep track of used genes
    gene_sites = {}  # Store sites for each gene

    with open(input_file, 'r') as file, \
         open(exceptions_file, 'w', newline='') as exceptions, \
         open(to_test_file, 'w', newline='') as to_test, \
         open(too_long_file, 'w', newline='') as too_long_file, \
         open(fasta_log_file, 'w', newline='') as fasta_log:

        reader = csv.reader(file, delimiter='\t')
        header = next(reader)  # Skip header

        exceptions_writer = csv.writer(exceptions, delimiter='\t')
        to_test_writer = csv.writer(to_test, delimiter='\t')
        too_long_writer = csv.writer(too_long_file, delimiter='\t')
        fasta_log_writer = csv.writer(fasta_log, delimiter='\t')

        exceptions_writer.writerow(header)
        to_test_writer.writerow(header)
        too_long_writer.writerow(header)
        fasta_log_writer.writerow(['Gene', 'Site(s)'])  # Write header for the FASTA log

        for row in reader:
            Kinase, gene, site, Score, alpha, beta, iso_B, alpha3, alphaI, hydturn, AFConf, NextAA, orga, site_aa, id_alphafold, in_vivo, in_vitro = row
            site = int(site)  # Convert site to an integer

            if gene not in gene_sites:
                gene_sites[gene] = []
            gene_sites[gene].append((site, site_aa, row))

        for gene, sites_data in gene_sites.items():
            sites_data = sorted(sites_data, key=lambda x: x[0])  # Sort by site number
            sites = [data[0] for data in sites_data]  # Extract the site numbers

            # Group sites based on proximity and process them accordingly
            processed_sites = set()  # Track processed sites
            for i, (site, site_aa, row) in enumerate(sites_data):
                if site in processed_sites:
                    continue  # Skip already processed sites
                
                # Group close sites (within 80 residues)
                close_sites = [site]
                for j in range(i + 1, len(sites_data)):
                    next_site, next_site_aa, next_row = sites_data[j]
                    if abs(site - next_site) <= 80:
                        close_sites.append(next_site)
                        processed_sites.add(next_site)
                    else:
                        break

                try:
                    # Fetch AlphaFold structure file
                    cif_file_path = os.path.join(cif_root_folder, f"AF-{row[14]}-F1-model_v4.cif.gz")
                    with gzip.open(cif_file_path, 'rt') as cif_file:
                        protein_sequence = extract_protein_sequence_from_cif(cif_file)
                        
                        if len(protein_sequence) > 800:  # Limit sequence length for FASTA output
                            too_long_writer.writerow(row)
                            protein_sequence = protein_sequence[max(0, site - 399): min(len(protein_sequence), site + 400)]

                    # Create phosphomimetic sequence
                    modified_sequence = phosphomimetic(protein_sequence, close_sites, site_aa)

                    # Write original and modified sequences to FASTA files
                    with open(os.path.join(output_folder, f"{gene}-E{site}.fasta"), 'w') as fasta_file:
                        fasta_file.write(f">{gene}-E{site}\n{modified_sequence}\n")
                    with open(os.path.join(output_folder, f"{gene}-{site}.fasta"), 'w') as fasta_file:
                        fasta_file.write(f">{gene}-{site}\n{protein_sequence}\n")

                    # Log the multi-site phosphorylation or the single site phosphorylation
                    if len(close_sites) > 1:
                        fasta_log_writer.writerow([gene, ','.join(map(str, close_sites))])
                    else:
                        fasta_log_writer.writerow([gene, str(site)])

                    to_test_writer.writerow(row)

                except Exception as e:
                    exceptions_writer.writerow(row)

    print(f"Processed {len(gene_sites)} genes.")

# Example usage
input_file = '/mnt/c/Users/crisd/Desktop/Src/Src_in_vitro_tested.txt'
cif_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/human_database"
output_folder = '/mnt/c/Users/crisd/Desktop/Src/fastas'
to_test_file = '/mnt/c/Users/crisd/Desktop/Src/fastas_created.txt'
exceptions_file = '/mnt/c/Users/crisd/Desktop/Src/exc.txt'
too_long_file = '/mnt/c/Users/crisd/Desktop/Src/toolong.txt'
fasta_log_file = '/mnt/c/Users/crisd/Desktop/Src/log.txt'
create_fastas(input_file, cif_folder, output_folder, to_test_file, exceptions_file, too_long_file, fasta_log_file)

print("FASTA file creation completed.")
