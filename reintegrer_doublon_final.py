import pandas as pd

# Charger les deux fichiers texte en DataFrames pandas
input_file = '/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_rmsd_filtered2.txt'
reference_file = '/mnt/c/Users/crisd/Downloads/Kinase_Substrate_Dataset/Kinase_Substrate_Dataset'

# Lire le fichier d'entrée unique
input_df = pd.read_csv(input_file, sep='\t')

# Lire le fichier de référence
reference_df = pd.read_csv(reference_file, sep='\t')

# Initialiser une liste pour stocker les nouvelles lignes
new_rows = []

# Parcourir chaque ligne du fichier d'entrée unique
for _, row in input_df.iterrows():
    gene = row['Gene']
    site = row['Res']
    
    # Chercher toutes les correspondances gene-site dans le fichier de référence
    matches = reference_df[(reference_df['SUB_GENE'] == gene) & (reference_df['SUB_MOD_RSD'].str.contains(str(site)))]
    
    for _, match in matches.iterrows():
        new_row = row.copy()
        new_row['kinase'] = match['KINASE']
        new_rows.append(new_row)

# Créer un nouveau DataFrame avec les lignes étendues
output_df = pd.DataFrame(new_rows)
output_df = output_df.sort_values(by='kinase')

# Sauvegarder le nouveau DataFrame dans un fichier texte
output_file = '/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_rmsd_filtered_final.txt'
output_df.to_csv(output_file, sep='\t', index=False)

print(f"Fichier étendu créé : {output_file}")
