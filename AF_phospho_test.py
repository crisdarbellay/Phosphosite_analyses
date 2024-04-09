import os
from main_tools2 import extract_residues_from_PDB
from main_tools2 import CompareTwoPDBs
import subprocess

def fatcat(cheminpdb1, cheminpdb2):
    # Commande FATCAT avec les chemins des fichiers PDB en arguments
    current_directory = os.getcwd()
    print("Répertoire actuel:", current_directory)          
    cmd = f"/home/crisdarbellay/stability_analyzer/FATCAT-dist/FATCATMain/FATCAT -b -p1 {cheminpdb1} -p2 {cheminpdb2} -o  output.aln"

    try:
        # Exécute la commande et récupère la sortie
        output = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
        
        # Divise la sortie en utilisant les espaces comme délimiteurs
        parts = output.stdout.split()
        
        # Recherche du RMSD par position dans la liste
        # Comme le RMSD est la 6ème valeur dans la liste
        rmsd = float(parts[6])
        return rmsd

    except:
        # En cas d'erreur lors de l'exécution de la commande
        print("Erreur lors de l'exécution de la commande FATCAT:")
        return None

# Dictionnaire de protéines et de leurs sites
protein_sites = {
        "glycogene_phosphorylase": "14",
        "BUB1":"969",
        "CDK2":"160",
        "ERK1":"185",
        "PAK1":"423",
        "RET":"905",
        "PKA":"197",
        "PKA":"338"
}

with open(r"/mnt/c/Users/crisd/Desktop/rmsd_comparisons.txt", "w") as f:
    for protein, site in protein_sites.items():
        os.chdir(f"/mnt/c/Users/crisd/Desktop/Protein_comparaison/{protein}")
        protein_folder = f"/mnt/c/Users/crisd/Desktop/Protein_comparaison/{protein}"

        if protein == "ERK1":
            real_unphospho_file = f"real-{protein}.pdb"
            real_phospho_file = f"real-phospho-ERK1.pdb"
            af_unphospho_file = f"AF-ERK1.pdb"
            af_d_phospho_file = f"ERK1-D183-D185.pdb"
            af_e_phospho_file = f"ERK1-E183-E185.pdb"

        else:
            real_unphospho_file = f"real-{protein}.pdb"
            real_phospho_file = f"real-phospho-{protein}-{site}.pdb"
            af_unphospho_file = f"{protein}-AF.pdb"
            af_d_phospho_file = f"{protein}-D{site}.pdb"
            af_e_phospho_file = f"{protein}-E{site}.pdb"
            sequence = extract_residues_from_PDB(real_unphospho_file)
            d = len(sequence)
            print(d)
            site = int(site)
            b = site - 1
            c = site + 1
            contigs = "A1-" + str(d)
        
        f.write(f"Protein {protein}, phosphosite {site}:\n")

        # AF unphosphorylated versus AF phosphomimetic structure D
        try:
                rmsd = CompareTwoPDBs(contigs, af_unphospho_file, af_d_phospho_file)
                f.write(f"\t AF unphosphorylated structure versus AF phosphomimetic structure D: {round(rmsd,2)}\n")
        except:
             pass

        # AF unphosphorylated versus AF phosphomimetic structure E
        try:
                rmsd = CompareTwoPDBs(contigs, af_unphospho_file, af_e_phospho_file)
                f.write(f"\t AF unphosphorylated structure versus AF phosphomimetic structure E: {round(rmsd,2)}\n")
        except:
             pass

        # Real unphosphorylated versus AF unphosphorylated
        try:
                rmsd = CompareTwoPDBs(contigs, af_unphospho_file, real_unphospho_file)
                f.write(f"\t real unphosphorylated structure versus AF unphosphorylated structure: {round(rmsd,2)}\n")
        except:
             pass

        # Real unphosphorylated versus AF phosphomimetic structure D
        try:
                rmsd = CompareTwoPDBs(contigs,af_d_phospho_file, real_unphospho_file)
                f.write(f"\t real unphosphorylated structure versus AF phosphomimetic structure D: {round(rmsd,2)}\n")
        except:
             pass

        # Real unphosphorylated versus AF phosphomimetic structure E
        try:
                rmsd = CompareTwoPDBs(contigs,af_e_phospho_file, real_unphospho_file)
                f.write(f"\t real unphosphorylated structure versus AF phosphomimetic structure E: {round(rmsd,2)}\n")
        except:
              pass

        # Real phosphorylated versus AF unphosphorylated
        try:
                rmsd = CompareTwoPDBs(contigs,af_unphospho_file, real_phospho_file)
                f.write(f"\t real phosphorylated structure versus AF unphosphorylated structure: {round(rmsd,2)}\n")
        except:
              pass

        # Real phosphorylated versus AF phosphomimetic structure D
        try:
                rmsd = CompareTwoPDBs(contigs,af_d_phospho_file, real_phospho_file)
                f.write(f"\t real phosphorylated structure versus AF phosphomimetic structure D: {round(rmsd,2)}\n")
        except:
              pass
        # Real phosphorylated versus AF phosphomimetic structure E
        try:
                rmsd = CompareTwoPDBs(contigs, af_e_phospho_file, real_phospho_file)
                f.write(f"\t real phosphorylated structure versus AF phosphomimetic structure E: {round(rmsd,2)}\n")
        except:
              pass
        # Real unphosphorylated versus real phosphorylated
        try:
                rmsd = CompareTwoPDBs(contigs, real_unphospho_file, real_phospho_file)
                f.write(f"\t real unphosphorylated structure versus real phosphorylated structure: {round(rmsd,2)}\n")
        except:
              pass