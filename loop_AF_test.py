import MDAnalysis as mda
import os
from main_tools import extract_residues_from_PDB

def calculate_rmsd(contigs,pdb1, pdb2):
    u1 = mda.Universe(pdb1, topology_format='PDB', format='PDB')
    u2 = mda.Universe(pdb2, topology_format='PDB', format='PDB')

    # Sélection des résidus à comparer en fonction des contigs
    selection = ""
    for contig in contigs.split("/"):
        chain, start, end = contig.split("-")
        selection += f"chain {chain} and resid {start}:{end} or "
    selection = selection.rstrip(" or ")

    ca1 = u1.select_atoms(f"name CA and ({selection})")
    ca2 = u2.select_atoms(f"name CA and ({selection})")

    # Calcul du RMSD
    rmsd = mda.analysis.rms.rmsd(ca1, ca2)
    return rmsd

def main():
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
            
            protein_folder = f"/mnt/c/Users/crisd/Desktop/Protein_comparaison/{protein}"

            if protein == "ERK1":
                real_unphospho_file = os.path.join(protein_folder, f"real-{protein}.pdb")
                real_phospho_file = os.path.join(protein_folder, f"real-phospho-ERK1.pdb")
                af_unphospho_file = os.path.join(protein_folder, f"AF-ERK1.pdb")
                af_d_phospho_file = os.path.join(protein_folder, f"ERK1-D183-D185.pdb")
                af_e_phospho_file = os.path.join(protein_folder, f"ERK1-E183-E185.pdb")
                contigs = "A1-"+str(b)+"/A"+str(c)+"-"+str(d)
            else:
                real_unphospho_file = os.path.join(protein_folder, f"real-{protein}.pdb")
                real_phospho_file = os.path.join(protein_folder, f"real-phospho-{protein}-{site}.pdb")
                af_unphospho_file = os.path.join(protein_folder, f"{protein}-AF.pdb")
                af_d_phospho_file = os.path.join(protein_folder, f"{protein}-D{site}.pdb")
                af_e_phospho_file = os.path.join(protein_folder, f"{protein}-E{site}.pdb")
                sequence = extract_residues_from_PDB(af_d_phospho_file)
                d = len(sequence)
                site = int(site)
                b = site - 1
                c = site + 1
                contigs = "A-1-"+str(d)

            print(contigs)
            f.write(f"Protein {protein}, phosphosite {site}:\n")

            # Real unphosphorylated versus real phosphorylated
            rmsd = calculate_rmsd(contigs, real_unphospho_file, real_phospho_file)
            f.write(f"\t real unphosphorylated structure versus real phosphorylated structure: {round(rmsd,2)}\n")

                # Real unphosphorylated versus AF unphosphorylated
            rmsd = calculate_rmsd(contigs, real_unphospho_file, af_unphospho_file)
            f.write(f"\t real unphosphorylated structure versus AF unphosphorylated structure: {round(rmsd,2)}\n")

                # Real unphosphorylated versus AF phosphomimetic structure D
            rmsd = calculate_rmsd(contigs, real_unphospho_file, af_d_phospho_file)
            f.write(f"\t real unphosphorylated structure versus AF phosphomimetic structure D: {round(rmsd,2)}\n")

                # Real unphosphorylated versus AF phosphomimetic structure E
            rmsd = calculate_rmsd(contigs, real_unphospho_file, af_e_phospho_file)
            f.write(f"\t real unphosphorylated structure versus AF phosphomimetic structure E: {round(rmsd,2)}\n")

                # Real phosphorylated versus AF unphosphorylated
            rmsd = calculate_rmsd(contigs, real_phospho_file, af_unphospho_file)
            f.write(f"\t real phosphorylated structure versus AF unphosphorylated structure: {round(rmsd,2)}\n")

                # Real phosphorylated versus AF phosphomimetic structure D
            rmsd = calculate_rmsd(contigs, real_phospho_file, af_d_phospho_file)
            f.write(f"\t real phosphorylated structure versus AF phosphomimetic structure D: {round(rmsd,2)}\n")

                # Real phosphorylated versus AF phosphomimetic structure E
            rmsd = calculate_rmsd(contigs, real_phospho_file, af_e_phospho_file)
            f.write(f"\t real phosphorylated structure versus AF phosphomimetic structure E: {round(rmsd,2)}\n")

                # AF unphosphorylated versus AF phosphomimetic structure D
            rmsd = calculate_rmsd(contigs, af_unphospho_file, af_d_phospho_file)
            f.write(f"\t AF unphosphorylated structure versus AF phosphomimetic structure D: {round(rmsd,2)}\n")

                # AF unphosphorylated versus AF phosphomimetic structure E
            rmsd = calculate_rmsd(contigs, af_unphospho_file, af_e_phospho_file)
            f.write(f"\t AF unphosphorylated structure versus AF phosphomimetic structure E: {round(rmsd,2)}\n")

            f.write("\n")

            
if __name__ == "__main__":
    main()
