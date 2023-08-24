import os

# Définir le dossier à parcourir
folder_path =  r"/mnt/c/Users/crisd/Desktop/ProteinDesign/results/mouse"

def remove_lines_with_none(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Filtrer les lignes ne contenant pas "None"
    filtered_lines = [line for line in lines if "100000" not in line and "None" not in line]

    with open(file_path, 'w') as file:
        file.writelines(filtered_lines)

# Parcourir les fichiers .txt dans le dossier
for file_name in os.listdir(folder_path):
    if file_name.endswith(".txt"):
        file_path = os.path.join(folder_path, file_name)
        remove_lines_with_none(file_path)

print("Suppression des lignes contenant 'None' terminée.")
