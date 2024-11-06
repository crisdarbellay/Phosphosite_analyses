import os

def remove_header_from_files(directory, header_text):
    # Parcourir tous les fichiers et sous-répertoires du répertoire
    for root, dirs, files in os.walk(directory):
        for filename in files:
            # Vérifier si c'est un fichier PDB
            if filename.endswith(".pdb"):
                file_path = os.path.join(root, filename)
                
                # Lire le contenu original du fichier
                with open(file_path, 'r') as file:
                    lines = file.readlines()
                
                # Filtrer les lignes pour supprimer celles qui correspondent au header
                new_lines = [line for line in lines if header_text not in line]
                
                # Écrire le nouveau contenu dans le fichier
                with open(file_path, 'w') as file:
                    file.writelines(new_lines)
                
                print(f"Lignes HEADER supprimées du fichier : {file_path}")

if __name__ == "__main__":
    # Le chemin du dossier qui contient les sous-dossiers avec les fichiers
    directory_path = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final_2"
    
    # Le texte exact du header à supprimer
    header_text = "HEADER                                            01-JUN-22"
    
    # Appeler la fonction pour supprimer les lignes avec le header
    remove_header_from_files(directory_path, header_text)
