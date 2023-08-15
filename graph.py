import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_data(input_file):
    with open(input_file, 'r') as fichier:
        lignes = fichier.readlines()

        en_tete = lignes[0].strip().split('\t')

        # Initialiser un dictionnaire de données
        tableau = {}

        # Parcourir les lignes à partir de la deuxième ligne
        for ligne in lignes[1:]:
            # Séparer chaque ligne en colonnes en utilisant le séparateur '\t'
            colonnes = ligne.strip().split('\t')
            
            # Extraire les informations clés
            gene = colonnes[0]
            site = int(colonnes[1])
            
            # Créer un dictionnaire de site s'il n'existe pas déjà
            if gene not in tableau:
                tableau[gene] = {}
            if site not in tableau[gene]:
                tableau[gene][site] = {}

            # Ajouter les données de chaque colonne au dictionnaire
            for i in range(2, len(colonnes)):
                colonne = en_tete[i]
                valeur = colonnes[i]
                tableau[gene][site][colonne] = valeur

        return tableau


def calculate_ranking(data, column_name, top_n=10, ascending=False):
    """
    Calculate ranking based on a specified column in the data dictionary.
    """

    # Create a list of tuples containing (gene, site, values) for the specified column
    ranking_data = []
    for gene, site_dict in data.items():
        for site, values in site_dict.items():
            if column_name in values:
                ranking_data.append((gene, site, values[column_name], values))

    # Sort the ranking data based on the specified column
    ranking_data.sort(key=lambda x: float(x[2]), reverse=not ascending)

    # Create a DataFrame for the ranking
    ranking_df = pd.DataFrame(ranking_data, columns=['Gene', 'Site', 'Confidence', 'Values'])
    ranking_df['Rank'] = range(1, 1 + len(ranking_df))

    return ranking_df.head(top_n)


def create_groups(data, input_file):
    with open(input_file, 'r') as fichier:
        first_line = fichier.readline()
        txt_columns = [col.strip() for col in first_line.split('\t') if col.endswith('.txt') or col.endswith('.txt\n')]
    
    groups = {
        'All Sites': [(gene, site) for gene, site_data in data.items() for site in site_data],
        'Not in any .txt': []
    }
    
    for gene, site_data in data.items():
        for site, values in site_data.items():
            found_in_txt = False
            for col in txt_columns:
                if values[col] == 'X':
                    found_in_txt = True
                    groups[col] = groups.get(col, [])
                    groups[col].append((gene, site))
                    break
            if not found_in_txt:
                groups['Not in any .txt'].append((gene, site))
        
    return groups


def calculate_group_stats(group_sites, data):
    confidence_total = 0
    score_total = 0
    count = len(group_sites)
    
    for site_info in group_sites:
        gene = site_info[0]
        site=site_info[1]
        site_data = data.get(gene, {}).get(site, {})
        confidence = float(site_data.get('Confidence', 0))
        score = float(site_data.get('Score', 0))
        confidence_total += confidence
        score_total += score
    
    if count == 0:
        return {
            'count': 0,
            'confidence_avg': 0.00,
            'score_avg': 0.00
        }
    
    confidence_avg = confidence_total / count
    score_avg = score_total / count
    
    return {
        'count': count,
        'confidence_avg': confidence_avg,
        'score_avg': score_avg
    }


def calculate_stats(data):
    stats = {}
    groups = data['Gene'].unique()

    for group in groups:
        group_data = data[data['Gene'] == group]
        count = len(group_data)
        confidence_avg = group_data['Confidence'].mean()
        score_avg = group_data['Score'].mean()
        
        stats[group] = {
            'count': count,
            'confidence_avg': confidence_avg,
            'score_avg': score_avg
        }
    
    return stats

def calculate_slices_stats(data, confidence_range):
    total_sites = 0
    phosphorylated_sites = 0
    total_score = 0
    total_nextAA = 0
    
    for gene_data in data.values():
        for site_data in gene_data.values():
            total_sites += 1  # Increment total_sites for each site
            
            confidence = float(site_data.get('Confidence', 0))
            if confidence_range[0] <= confidence < confidence_range[1]:
                phosphorylated_sites += 1
                total_score += float(site_data.get('Score', 0))
                total_nextAA += float(site_data.get('NextAA', 0))
    
    if total_sites == 0:
        return {
            'Number of phosphorylated site': 0,
            'Frequency': '0.00%',
            'Average score': 0.00,
            'Average nextAA in sec struct': 0.00
        }
    
    frequency = (phosphorylated_sites / total_sites) * 100
    if not phosphorylated_sites==0:
        avg_score = total_score / phosphorylated_sites
        avg_nextAA = total_nextAA / phosphorylated_sites
    else:
        avg_nextAA=None
        avg_score=None
    
    return {
        'Number of phosphorylated site': phosphorylated_sites,
        'Frequency': f'{frequency:.2f}%',
        'Average score': avg_score,
        'Average nextAA in sec struct': avg_nextAA
    }

def calculate_all_slices_stats(data):
    confidence_ranges = [(20, 30), (30, 40), (40, 50), (50, 60), (60, 70), (80, 90)]
    slices_stats = {}
    
    for confidence_range in confidence_ranges:
        slice_name = f'Slice of confidence score: {confidence_range[0]}-{confidence_range[1]}'
        slices_stats[slice_name] = calculate_slices_stats(data, confidence_range)
    
    return slices_stats

def create_heatmap(data, group_sites, group_name, characteristic, output_folder):
    plt.figure(figsize=(10, 6))
    
    max_confidence = 100  # La plage des valeurs de confiance est fixe de 0 à 100
    max_value = 0  # Initialiser à 0

    for gene, site in group_sites:
        site_data = data.get(gene, {}).get(site, {})
        confidence = float(site_data.get('Confidence', 0))
        value = float(site_data.get(characteristic, 0))

        max_value = max(max_value, value)  # Mettre à jour la valeur maximale

    heatmap_data = np.zeros((101, int(max_value) + 1))

    for gene, site in group_sites:
        site_data = data.get(gene, {}).get(site, {})
        confidence = float(site_data.get('Confidence', 0))
        value = float(site_data.get(characteristic, 0))

        heatmap_data[int(confidence)][int(value)] += 1

    heatmap_data_transposed = np.transpose(heatmap_data)  # Transposer la matrice

    plt.imshow(heatmap_data_transposed, cmap='coolwarm', origin='lower', extent=[0, max_confidence, 0, max_value], aspect='auto', interpolation='bilinear')
    plt.colorbar(label='Frequency')
    plt.xlabel('Alphafold Confidence')
    plt.ylabel(characteristic)
    plt.title(f'Heatmap of {characteristic} vs. Alphafold Confidence - Group: {group_name}')
    plt.savefig(f'{output_folder}/{group_name}_{characteristic}_heatmap.png')
    plt.show()




def print_outside_file(output_file_path,input_file,to_rank,secondary_structures):

    data = parse_data(input_file)
    groups = create_groups(data, input_file)
    datas=calculate_group_stats(groups, data)
    slices=calculate_all_slices_stats(data)

    with open(input_file) as fichier:
        first_line = fichier.readline()
        txt_columns = [col.strip() for col in first_line.split('\t') if col.endswith('.txt') or col.endswith('.txt\n')]

    with open(output_file_path, 'w') as output_file:

        for group_name, group_sites in groups.items():
            group_stats = calculate_group_stats(group_sites, data)
            output_file.write(f"Group {group_name}\n")
            output_file.write(f"Number of sites : {group_stats['count']}\nAverage alphafold confidence : {group_stats['confidence_avg']:.2f}\n")
            output_file.write(f"Average score : {group_stats['score_avg']:.2f}\n\n")
            for characteristic in ['Score','NextAA']:
                create_heatmap(data, group_sites, group_name, characteristic, output_folder)

        for rank in to_rank:
            ranking = calculate_ranking(data, rank[0], rank[1], rank[2])
            output_file.write(f"Ranking {rank[0]}:\n")
            output_file.write("\tGene\tsite\tScore\tAF_Conf\tNextAA\t")
            for secondary_structure in secondary_structures:
                output_file.write(f"{secondary_structure[2]}\t")
            output_file.write("\n")
            for index, row in ranking.iterrows():
                output_file.write(f"{index + 1}.\t{row['Gene']}\t{row['Site']}\t{row['Values']['Score']}\t{row['Values']['Confidence']}\t{row['Values']['NextAA']}\t")
                for secondary_structure in secondary_structures:
                    output_file.write(f"{row['Values'].get(secondary_structure[2], '')}\t")
                output_file.write("\n")
            output_file.write("\n")

        for slice_name, slice_data in slices.items():
            output_file.write(f"{slice_name}:\n")
            for key, value in slice_data.items():
                output_file.write(f"{key}: {value}\n")
            output_file.write("\n")

input_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/PKA/results.txt" 
output_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/PKA/results"

alpha_H_distance=10
beta_E_distance=10
count=0
secondary_structures=[('H',alpha_H_distance,'alpha',count),('E',beta_E_distance,'beta',count),('B',beta_E_distance,'iso_B',count),('G',alpha_H_distance,'alpha3',count),('I',alpha_H_distance,'alphaI',count),('T',alpha_H_distance,'hydrogene_turn',count)]
to_rank=[('Score',10,False),('Confidence',10,False),('NextAA',10,True),('alpha',10,False),('beta',10,False),('iso_B',10,False),('alpha3',10,False),('alphaI',10,False)]
                
output_file_path = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/PKA/results/output_file.txt"  # Replace with the desired output file path

#print_outside_file(output_file_path,input_file,to_rank,secondary_structures)