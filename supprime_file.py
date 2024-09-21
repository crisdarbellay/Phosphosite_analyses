import os

# Chemin du dossier principal
dossier_principal = '/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_1/pdb-ball-30-contigs-adapted'

# Parcours récursif de tous les sous-dossiers
for root, dirs, files in os.walk(dossier_principal):
    for file in files:
        if file.startswith('no_oxt'):
            chemin_fichier = os.path.join(root, file)
            os.remove(chemin_fichier)
            print(f"Supprimé: {chemin_fichier}")

print("Opération terminée.")
