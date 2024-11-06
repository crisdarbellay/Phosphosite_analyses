import os

def add_header_to_files(directory, header_text):
    # Parcourir tous les fichiers et sous-répertoires du répertoire
    for root, dirs, files in os.walk(directory):
        for filename in files:
            # Vérifier si c'est un fichier PDB
            if filename.endswith(".pdb"):
                file_path = os.path.join(root, filename)
                
                # Lire le contenu original du fichier
                with open(file_path, 'r') as file:
                    original_content = file.read()
                    if original_content.split('\n')[0] == header_text:
                        print('header deja présent pour',file)
                        continue
                
                # Créer le nouveau contenu avec le header
                new_content = header_text + "\n" + original_content
                
                # Écrire le nouveau contenu dans le fichier
                with open(file_path, 'w') as file:
                    file.write(new_content)
                
                print(f"Header ajouté au fichier : {file_path}")

if __name__ == "__main__":
    # Le chemin du dossier qui contient les sous-dossiers avec les fichiers
    directory_path = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final_2"
    
    # Le texte du header exact à ajouter en haut de chaque fichier
    header_text = "HEADER                                            01-JUN-22"
    
    # Appeler la fonction pour ajouter le header à chaque fichier PDB
    add_header_to_files(directory_path, header_text)
