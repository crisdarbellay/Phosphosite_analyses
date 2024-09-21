import os
import random
import csv
import gzip

#This code generate control to test the effect of the phosphomimetic and to test the selected proteins

def read_input_file(input_file):
    with open(input_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)
        data = [row for row in reader]
    return header, data

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

def mutate_sequence(sequence, site, new_aa):
    sequence_list = list(sequence)
    sequence_list[site - 1] = new_aa
    return ''.join(sequence_list)

def generate_fasta_files(data, cif_root_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    log_entries = []

    for row in data:
        kinase, gene, site, score, percent, alpha, beta, iso_B, alpha3, alphaI, hydturn, AFConf, NextAA, orga, site_aa, sub_id, in_vivo, in_vitro = row
        site = int(site)

        try:
            cif_folder = os.path.join(cif_root_folder, "human_database")
            cif_file_path = os.path.join(cif_folder, f"AF-{sub_id}-F1-model_v4.cif.gz")
            with gzip.open(cif_file_path, 'rt') as cif_file:
                protein_sequence = extract_protein_sequence_from_cif(cif_file)
                len_prot = len(protein_sequence)
                if len_prot > 800:
                    site_idx = site - 1
                    start = max(site_idx - 399, 0)
                    end = min(site_idx + 400, len(protein_sequence))
                    protein_sequence = protein_sequence[start:end]
                    fasta_filename = f"{gene}-{site}"
                    if site_idx < 400:
                        site = site
                    else:
                        site = 400
                else:
                    fasta_filename = gene
        except Exception as e:
            continue

        # Control 1: Phosphosite mutated to A
        mutated_sequence_1 = mutate_sequence(protein_sequence, site, 'A')
        with open(os.path.join(output_folder, f"{fasta_filename}_control1.fasta"), 'w') as fasta_file:
            fasta_file.write(f">{fasta_filename}_control1\n{mutated_sequence_1}\n")
        log_entries.append([kinase, gene, site, site, 'S/T/Y', 'A', score, percent, alpha, beta, iso_B, alpha3, alphaI, hydturn, AFConf, NextAA, orga, site_aa, sub_id, in_vivo, in_vitro])

        # Control 2: Random S mutated to D, another to E
        s_sites = [i for i, aa in enumerate(protein_sequence) if aa == 'S' and i != site - 1]
        if len(s_sites) >= 2:
            random_sites = random.sample(s_sites, 2)
            mutated_sequence_2_d = mutate_sequence(protein_sequence, random_sites[0] + 1, 'D')
            mutated_sequence_2_e = mutate_sequence(mutated_sequence_2_d, random_sites[1] + 1, 'E')
            with open(os.path.join(output_folder, f"{fasta_filename}_control2.fasta"), 'w') as fasta_file:
                fasta_file.write(f">{fasta_filename}_control2\n{mutated_sequence_2_e}\n")
            log_entries.append([kinase, gene, random_sites[0] + 1, random_sites[0] + 1, 'S', 'D', score, percent, alpha, beta, iso_B, alpha3, alphaI, hydturn, AFConf, NextAA, orga, site_aa, sub_id, in_vivo, in_vitro])
            log_entries.append([kinase, gene, random_sites[1] + 1, random_sites[1] + 1, 'S', 'E', score, percent, alpha, beta, iso_B, alpha3, alphaI, hydturn, AFConf, NextAA, orga, site_aa, sub_id, in_vivo, in_vitro])

        # Control 3: Random AA mutated to D, another to E
        all_sites = [i for i, aa in enumerate(protein_sequence) if i != site - 1]
        if len(all_sites) >= 2:
            random_sites = random.sample(all_sites, 2)
            mutated_sequence_3_d = mutate_sequence(protein_sequence, random_sites[0] + 1, 'D')
            mutated_sequence_3_e = mutate_sequence(mutated_sequence_3_d, random_sites[1] + 1, 'E')
            with open(os.path.join(output_folder, f"{fasta_filename}_control3.fasta"), 'w') as fasta_file:
                fasta_file.write(f">{fasta_filename}_control3\n{mutated_sequence_3_e}\n")
            log_entries.append([kinase, gene, random_sites[0] + 1, random_sites[0] + 1, protein_sequence[random_sites[0]], 'D', score, percent, alpha, beta, iso_B, alpha3, alphaI, hydturn, AFConf, NextAA, orga, site_aa, sub_id, in_vivo, in_vitro])
            log_entries.append([kinase, gene, random_sites[1] + 1, random_sites[1] + 1, protein_sequence[random_sites[1]], 'E', score, percent, alpha, beta, iso_B, alpha3, alphaI, hydturn, AFConf, NextAA, orga, site_aa, sub_id, in_vivo, in_vitro])

    return log_entries

def write_log_file(log_entries, log_file_path):
    with open(log_file_path, 'w', newline='') as log_file:
        writer = csv.writer(log_file, delimiter='\t')
        writer.writerow(["kinase", "gene", "contrôlé", "testé", "origine", "nouveau", "Score", "percent", "alpha", "beta", "iso_B", "alpha3", "alphaI", "hydturn", "AFConf", "NextAA", "orga", "site_aa", "sub_id", "in_vivo", "in_vitro"])
        for entry in log_entries:
            writer.writerow(entry)

# Utilisation
input_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/too_long.txt"
cif_root_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign"
output_folder = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/controls/fastas"
log_file_path = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/controls/log_file.txt"

header, data = read_input_file(input_file)
log_entries = generate_fasta_files(data, cif_root_folder, output_folder)
write_log_file(log_entries, log_file_path)

print("Contrôles générés avec succès.")
