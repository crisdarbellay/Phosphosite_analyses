from collections import OrderedDict

selection_path = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/to_test/selection_full.txt"

# Lire les lignes du fichier et enlever les doublons tout en conservant l'ordre
with open(selection_path, 'r') as file:
    lines = file.readlines()
    unique_lines = list(OrderedDict.fromkeys(lines))

# Écrire les lignes uniques dans le fichier
with open(selection_path, 'w') as file:
    for line in unique_lines:
        file.write(line)
