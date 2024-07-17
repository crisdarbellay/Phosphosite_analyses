import os
import re

def rename_protein_files(directory):
    # Dictionnaire pour stocker les associations protein -> site
    protein_site_map = {}

    # Parcours des fichiers dans le dossier
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            match_d = re.match(r"^(.*?)-D(\d+)\.fasta$", filename)
            if match_d:
                protein = match_d.group(1)
                site = int(match_d.group(2))
                if protein in protein_site_map:
                    protein_site_map[protein] = min(protein_site_map[protein], site)
                else:
                    protein_site_map[protein] = site

    # Renommer les fichiers {protein}.fasta en {protein}-{site}.fasta
    for filename in os.listdir(directory):
        if filename.endswith(".fasta") and not re.match(r"^(.*?)-[DE]\d+\.fasta$", filename):
            protein = filename[:-6]  # Enlever l'extension .fasta
            if protein in protein_site_map:
                new_filename = f"{protein}-{protein_site_map[protein]}.fasta"
                os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))
                print(f"Renamed {filename} to {new_filename}")
            else:
                print(f"No site found for {filename}, skipping.")

# Exemple d'utilisation
rename_protein_files("/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/fastas-ball-30")
