def filter_results(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extraire le header
    header = lines[0]

    # Liste pour stocker les résultats filtrés
    filtered_lines = [header]

    # Parcourir chaque ligne de données
    for line in lines[1:]:
        # Séparer les valeurs par tabulation
        columns = line.split('\t')

        # Extraire les valeurs des colonnes 18, 19, 20, 21
        values = [float(columns[18]), float(columns[19]), float(columns[20]), float(columns[21])]

        # Vérifier si une des valeurs est supérieure à 4
        if any(value > 4 for value in values):
            filtered_lines.append(line)

    return filtered_lines

def save_filtered_results(filtered_lines, output_file):
    with open(output_file, 'w') as file:
        for line in filtered_lines:
            file.write(line)

# Utilisation des fonctions
input_file = '/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/3_algorithm/results/results_rmsd.txt'
output_file = '/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/3_algorithm/results/filtered_results.txt'

filtered_lines = filter_results(input_file)
save_filtered_results(filtered_lines, output_file)

print(f"Filtered results saved to {output_file}")
