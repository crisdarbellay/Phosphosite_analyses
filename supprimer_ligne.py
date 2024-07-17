selection_path = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/to_test/selection_full.txt"

def line_condition(line):
    values = line.strip().split()
    if len(values) < 2:
        # Ignore the line if it does not have enough columns
        return False
    try:
        # Vérifie si la dernière ou avant-dernière valeur est supérieure à 2
        last_value = float(values[-1])
        second_last_value = float(values[-2])
        return last_value > 2 or second_last_value > 2
    except ValueError:
        # Si la conversion en float échoue (par exemple pour les en-têtes), ignorer la ligne
        return False

# Lire les lignes du fichier et enlever celles ne satisfaisant pas la condition
with open(selection_path, 'r') as file:
    lines = file.readlines()

# Sauvegarder l'en-tête
header = lines[0]
filtered_lines = [header] + [line for line in lines[1:] if line_condition(line)]

# Écrire les lignes filtrées dans le fichier
with open(selection_path, 'w') as file:
    for line in filtered_lines:
        file.write(line)
