import pandas as pd
from prettytable import PrettyTable

# Chemin du fichier d'entrée
input_file = '/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_rmsd_filtered_final.txt'  # Remplacez par le chemin de votre fichier d'entrée

# Chemin du fichier de sortie
output_file = '/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_rmsd_filtered_in_vitro_pretty.txt'  # Remplacez par le chemin de votre fichier de sortie

# Lire le fichier d'entrée en DataFrame pandas
df = pd.read_csv(input_file, sep='\t')

# Filtrer les lignes où la colonne 'in_vitro' contient 'X'
filtered_df = df[df['in_vitro'] == 'X']

# Extraire l'en-tête et les données filtrées
header = filtered_df.columns.tolist()
results_filtered = filtered_df.values.tolist()

# Fonction pour sauvegarder les résultats triés dans un fichier avec PrettyTable
def write_pretty_results(file_path, header, results_filtered, use_pretty_table=True):
    if use_pretty_table:
        table = PrettyTable()
        table.field_names = header
        for result in results_filtered:
            table.add_row(result)
        table.border = False
        table.header = True
        with open(file_path, 'w') as f_out:
            f_out.write(table.get_string())
    else:
        with open(file_path, 'w') as f_out:
            f_out.write('\t'.join(header) + '\n')
            for result in results_filtered:
                f_out.write('\t'.join(map(str, result)) + '\n')

# Sauvegarder le DataFrame filtré dans un fichier texte joli
write_pretty_results(output_file, header, results_filtered)

print(f"Fichier filtré et joli créé : {output_file}")
