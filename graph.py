import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def create_heatmap(x, y, z, title, output_path):
    plt.figure(figsize=(10, 6))
    cmap = cm.get_cmap('coolwarm')
    plt.imshow(z, cmap=cmap, extent=[x.min(), x.max(), y.min(), y.max()], origin='lower', aspect='auto')
    
    plt.colorbar(label='Frequency (%)')
    plt.xlabel('Score')
    plt.ylabel('Alphafold confidence score')
    plt.title(title)
    
    plt.savefig(output_path)
    plt.close()

def get_top_n_residues(group_data, n, sort_key_index, top=False):
    # Triez la liste group_data en fonction du critère de classement (sort_key_index)
    if top:
        sorted_data = sorted(group_data, key=lambda data: data[sort_key_index], reverse=True)[::-1]
    else:    
        sorted_data = sorted(group_data, key=lambda data: data[sort_key_index], reverse=True)
    
    # Sélectionnez les N premiers éléments (top N) de la liste triée
    top_n_residues = sorted_data[:n]
    
    return top_n_residues

def create_surface_plot(x, y, z, title, output_path):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    
    ax.set_xlabel('Score')
    ax.set_ylabel('Alphafold confidence score')
    ax.set_zlabel('Frequency (%)')
    ax.set_title(title)
    
    fig.colorbar(surf, shrink=0.5, aspect=10)
    
    plt.savefig(output_path)
    plt.close()

def write_ranking_data(output_file, title, data_list):
    output_file.write(f"{title}:\n\tGene\tsite\tScore\tAlphafold confidence\tNextAA in sec struct\n")
    for i, data in enumerate(data_list, start=1):
        confidence, score, next_aa, gene, site = data
        output_file.write(f"{i}.\t{gene}\t{site}\t{score:.2f}\t\t{confidence:.2f}\t\t{next_aa:.0f}\n")
    output_file.write("\n")

def write_slice_info(output_file, confidence_min, confidence_max, total_sites, filtered_data):
    output_file.write(f"Slice of confidence score: {confidence_min}-{confidence_max}\n")
    output_file.write(f"Number of phosphorylated site: {len(filtered_data)}\n")
    output_file.write(f"Frequency: {len(filtered_data) / total_sites * 100:.2f}%\n")
    output_file.write(f"Average score: {np.mean([data[1] for data in filtered_data]):.2f}\n")
    output_file.write(f"Average nextAA in sec struct: {np.mean([data[2] for data in filtered_data]):.2f}\n\n")

def generate_3d_surface_graphs(input_file, output_folder):
    # Créer un dossier de destination s'il n'existe pas
    os.makedirs(output_folder, exist_ok=True)
    
    # Lire le fichier ligne par ligne et traiter les données
    with open(input_file, 'r') as file:
        # Ignorer la première ligne (en-tête)
        next(file)
        
        # Initialiser des dictionnaires pour chaque groupe
        groups = {
            'albuquerque_data': [],
            'soulard_data': [],
            'swaney_data': [],
            'in_no_data': [],
            'in_all_data': [],
            'all_residue': []  # Ajout du groupe "all_residue"
        }
        
        for line in file:
            parts = line.split('\t')
            gene = parts[0]
            site = parts[1]
            score = float(parts[2])
            if parts[5] != '':
                confidence = float(parts[5])
                next_aa = float(parts[7])
            alb_txt = parts[8]
            sou_txt = parts[9]
            swa_txt = parts[10]
            
            # Ajouter les données au groupe approprié en fonction des colonnes Alb.txt, Sou.txt et Swa.txt
            if alb_txt == 'X':
                groups['albuquerque_data'].append((confidence, score, next_aa,gene,site))
            if sou_txt == 'X':
                groups['soulard_data'].append((confidence, score, next_aa,gene,site))
            if swa_txt == 'X':
                groups['swaney_data'].append((confidence, score, next_aa,gene,site))
            if alb_txt != 'X' and sou_txt != 'X' and swa_txt == 'X':
                groups['in_no_data'].append((confidence, score, next_aa,gene,site))
            if alb_txt == 'X' and sou_txt == 'X' and swa_txt == 'X':
                groups['in_all_data'].append((confidence, score, next_aa,gene,site))
            
            # Ajouter les données au groupe "all_residue"
            groups['all_residue'].append((confidence, score, next_aa,gene,site))
    
    # Générer les graphiques en forme de surface pour chaque groupe
    for group_name, group_data in groups.items():
        if not group_data:
            continue
        
        # Créer les tranches de score de confiance
        confidence_bins = np.arange(0, 101, 2)
        score_bins = np.arange(0, 23, 1)  # Crée des tranches de 2 en 2 jusqu'à 22
        confidence_bins10 = np.arange(0, 101, 10)

        # Compter le nombre d'occurrences dans chaque tranche
        freq_matrix = np.zeros((len(confidence_bins) - 1, len(score_bins) - 1))
        
        for data in group_data:
            confidence = data[0]
            score = data[1]
            
            confidence_bin = np.digitize(confidence, confidence_bins) - 1
            score_bin = np.digitize(score, score_bins) - 1
            
            if confidence_bin < 0 or confidence_bin >= len(confidence_bins) - 1 or score_bin < 0 or score_bin >= len(score_bins) - 1:
                continue
            
            freq_matrix[confidence_bin, score_bin] += 1
        
        # Convertir les occurrences en fréquences en pourcentage
        total_sites = len(group_data)
        freq_matrix = freq_matrix / total_sites * 100
        
        create_heatmap(
            *np.meshgrid(score_bins[:-1], confidence_bins[:-1]),
            freq_matrix,
            f'{group_name.capitalize()} - Frequency in function of Alphafold confidence and number of AA in a ball within 10 angstrom',
            os.path.join(output_folder, f'{group_name}_heatmap.png')
        )
        
    # Écrire les informations dans le fichier de sortie
    output_filename = os.path.join(output_folder, f'{group_name}_info.txt')
    with open(output_filename, 'w') as output_file:
        ecrire_infos_sortie(output_file, groups, confidence_bins10, total_sites)
    print(f"Informations écrites")

def ecrire_infos_sortie(output_file, groupes, tranches_confiance, total_sites):
    for nom_groupe, donnees_groupe in groupes.items():
        if not donnees_groupe:
            continue

        output_file.write(f"Groupe {nom_groupe}\n")
        output_file.write(f"Nombre de résidus : {len(donnees_groupe)}\t")
        output_file.write(f"Moyenne de confiance Alphafold : {np.mean([data[0] for data in donnees_groupe]):.2f}\n")
        output_file.write(f"Moyenne de score : {np.mean([data[1] for data in donnees_groupe]):.2f}\n\n")
        
    top_10_score_residues = get_top_n_residues(groupes['all_residue'], 10, sort_key_index=1)
    top_10_confidence_residues = get_top_n_residues(groupes['all_residue'], 10, sort_key_index=0)
    top_10_nextaa_residues = get_top_n_residues(groupes['all_residue'],10,sort_key_index=2,top=True)

    write_ranking_data(output_file, "Ranking score", top_10_score_residues)
    write_ranking_data(output_file, "Ranking Alphafold confidence", top_10_confidence_residues)
    write_ranking_data(output_file, "Ranking NextAA in sec struct", top_10_nextaa_residues)
    
    for confidence_bin in range(len(tranches_confiance) - 1):
        confidence_min = tranches_confiance[confidence_bin]
        confidence_max = tranches_confiance[confidence_bin + 1]
        
        filtered_data = [data for data in groupes['all_residue'] if confidence_min <= data[0] <= confidence_max]
        if filtered_data:
            write_slice_info(output_file, confidence_min, confidence_max, total_sites, filtered_data)

# Appeler la fonction pour générer les graphiques en forme de surface
input_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/results.txt"  # Remplacez par le chemin de votre fichier d'entrée
output_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/Cdk1/Graph"
generate_3d_surface_graphs(input_file, output_folder)
