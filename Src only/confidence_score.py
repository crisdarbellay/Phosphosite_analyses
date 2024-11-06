import os
from Bio.PDB import PDBParser
import math
import csv

def calculate_distance(coord1, coord2):
    """Calcule la distance Euclidienne entre deux coordonnées 3D."""
    return math.sqrt(sum((c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2)))

def find_pdb_file(gene, site, variant, base_dir):
    """Trouve le fichier PDB dans le sous-dossier correspondant au gène et au site."""
    if variant == "":
        # Chercher soit {gene}-{site} soit {gene} pour le wild type
        subfolders = [f"{gene}-{site}", f"{gene}"]
    else:
        # Chercher les variants {gene}-D{site} ou {gene}-E{site}
        subfolders = [f"{gene}-{variant}{site}"]

    for subfolder in subfolders:
        pdb_dir = os.path.join(base_dir, subfolder)
        if os.path.exists(pdb_dir):
            for pdb_file in os.listdir(pdb_dir):
                if pdb_file.endswith(".pdb"):
                    return os.path.join(pdb_dir, pdb_file)
    return None

def extract_confidence_score(pdb_file, site_of_interest):
    """Extrait le score de confiance pour le résidu à proximité d'un site d'intérêt."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]  # Assume qu'il n'y a qu'un modèle dans la structure

    # Trouver le résidu d'intérêt (celui du site donné)
    target_residue = None
    target_coord = None
    for chain in model:
        for residue in chain.get_residues():
            if residue.get_id()[1] == site_of_interest:
                target_residue = residue
                target_coord = residue['CA'].get_coord()  # Coordonnées de l'atome CA
                break
        if target_residue:
            break
    
    if not target_residue:
        print(f"Résidu {site_of_interest} introuvable dans {pdb_file}")
        return None
    
    # Calculer les distances CA-CA avec les autres résidus et extraire les scores de confiance
    confidence_scores = []
    for chain in model:
        for residue in chain.get_residues():
            if 'CA' in residue:
                coord = residue['CA'].get_coord()
                distance = calculate_distance(target_coord, coord)
                if distance <= 6.0:  # Seuil de 6 Å
                    confidence_score = residue['CA'].bfactor  # Le score de confiance est dans le bfactor
                    confidence_scores.append(confidence_score)
    
    # Calculer et retourner la moyenne des scores de confiance
    if confidence_scores:
        return sum(confidence_scores) / len(confidence_scores)
    return None

def process_input_file(input_file, base_pdb_dir, output_file):
    """Parcourt le fichier input, calcule les différences de confiance et écrit les résultats."""
    with open(input_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)  # Lire l'en-tête
        rows = list(reader)  # Lire toutes les lignes

    # Ajouter trois nouvelles colonnes à l'en-tête
    header.extend(['Diff_WT-D', 'Diff_WT-E', 'Diff_D-E'])

    results = []

    # Parcourir les lignes du fichier
    for row in rows:
        if row[0] !='Src' and row[0]!='PKACA':
            continue
        gene = row[1]  # Colonne 1: Gene
        site = int(row[2])  # Colonne 2: Site (numéro)

        # Rechercher le fichier PDB pour le wild type
        try:
            pdb_file = find_pdb_file(gene, site, "", base_pdb_dir)
            confidence_score_wt = extract_confidence_score(pdb_file, site) if pdb_file else None
        except Exception as e:
            print(f"Erreur lors de la recherche du wild type pour {gene}-{site}: {e}")
            confidence_score_wt = None

        # Rechercher les fichiers PDB pour les variants D et E
        pdb_file_D = find_pdb_file(gene, site, "D", base_pdb_dir)
        confidence_score_D = extract_confidence_score(pdb_file_D, site) if pdb_file_D else None

        pdb_file_E = find_pdb_file(gene, site, "E", base_pdb_dir)
        confidence_score_E = extract_confidence_score(pdb_file_E, site) if pdb_file_E else None

        # Calculer les différences de confiance
        diff_wt_d = round(confidence_score_wt - confidence_score_D,1) if confidence_score_wt and confidence_score_D else None
        diff_wt_e = round(confidence_score_wt - confidence_score_E,1) if confidence_score_wt and confidence_score_E else None
        diff_d_e = round(confidence_score_D - confidence_score_E,1) if confidence_score_D and confidence_score_E else None

        # Ajouter les différences aux lignes
        row.extend([diff_wt_d, diff_wt_e, diff_d_e])
        results.append(row)

    # Écrire les résultats dans un nouveau fichier avec les colonnes supplémentaires
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(header)  # Écrire l'en-tête
        writer.writerows(results)  # Écrire les lignes avec les différences de confiance

# Exécuter le traitement
input_file = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/8_angstrom/2_algorithm/phosphosites/results/results_rmsd_score.txt"  # Spécifie le chemin du fichier d'entrée
base_pdb_dir = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/8_angstrom/2_algorithm/mass_colab_full/pdb"  # Dossier de base contenant les sous-dossiers PDB
output_file = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/8_angstrom/2_algorithm/phosphosites/results/results_with_confidence_diff.txt"  # Chemin du fichier de sortie

# process_input_file(input_file, base_pdb_dir, output_file)
