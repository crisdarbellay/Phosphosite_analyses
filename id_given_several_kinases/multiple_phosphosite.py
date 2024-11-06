import pandas as pd
import os

def extract_multiple_in_vitro_phosphorylations(file_path, output_directory):
    # Lire le fichier CSV
    df = pd.read_csv(file_path, sep='\t', encoding='iso-8859-1')

    # Filtrer pour Src humain et phosphorylations in vitro uniquement
    filtered_df = df[(df['KINASE'] == 'Src') & 
                     (df['KIN_ORGANISM'] == 'human') & 
                     (df['IN_VITRO_RXN'] == 'X')]

    # Compter les phosphorylations *in vitro* par gène
    phosphorylation_counts = filtered_df.groupby('SUB_GENE').size()

    # Sélectionner les gènes phosphorylés plusieurs fois *in vitro*
    multiple_phosphorylated_genes = phosphorylation_counts[phosphorylation_counts > 1].index

    # Filtrer à nouveau pour ne garder que ces gènes
    final_df = filtered_df[filtered_df['SUB_GENE'].isin(multiple_phosphorylated_genes)]

    # Regrouper par gène et concaténer uniquement les sites *in vitro* (SUB_MOD_RSD)
    grouped_df = final_df.groupby('SUB_GENE').agg({
        'KINASE': 'first',
        'SUB_MOD_RSD': lambda x: ' '.join(x),
        'SITE_+/-7_AA': lambda x: ' | '.join(x),  # Si tu veux aussi afficher les séquences +/-7 AA, sinon supprime cette ligne
        'SUB_ACC_ID': 'first'
    }).reset_index()

    # Créer un fichier de sortie
    output_file = os.path.join(output_directory, "src_multiple_in_vitro_phosphorylated_genes_combined.txt")

    # Sauvegarder les résultats avec seulement les phosphorylations *in vitro*
    grouped_df[['KINASE', 'SUB_GENE', 'SUB_MOD_RSD', 'SUB_ACC_ID']].to_csv(output_file, sep='\t', index=False)

    print(f"Results saved to {output_file}")

# Exemple d'utilisation :
file_path = r"/mnt/c/Users/crisd/Downloads/Kinase_Substrate_Dataset/Kinase_Substrate_Dataset"
output_directory = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/4_algorithm/test_final"

extract_multiple_in_vitro_phosphorylations(file_path, output_directory)
