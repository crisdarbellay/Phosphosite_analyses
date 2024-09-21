import pandas as pd

def filter_organism(input_file, output_file):
    # Lire le fichier d'entrée
    df = pd.read_csv(input_file, delimiter='\t')

    # Filtrer les lignes où l'organisme n'est pas humain, rat ou souris
    valid_organisms = ['human', 'rat', 'mouse']
    df_filtered = df[df['KIN_ORGANISM'].isin(valid_organisms) & df['SUB_ORGANISM'].isin(valid_organisms)]

    # Écrire les résultats dans un nouveau fichier
    df_filtered.to_csv(output_file, sep='\t', index=False)

# Exemple d'utilisation

output_file = r"/mnt/c/Users/crisd/Downloads/Kinase_Substrate_Dataset/Kinase_Substrate_Dataset_final"
input_file = r"/mnt/c/Users/crisd/Downloads/Kinase_Substrate_Dataset/Kinase_Substrate_Dataset_unique_no_doublon"

filter_organism(input_file, output_file)

print("Processing completed.")

