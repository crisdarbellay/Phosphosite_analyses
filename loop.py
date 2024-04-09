from Bio.PDB import Superimposer
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import os
from main_tools import extract_residues_from_PDB

three_to_one_letter = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V'
}

def extract_gene_site_pairs(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    gene_site_pairs = []

    section_keywords = ["Ranking Score", "Ranking AFConf", "Ranking NextAA"]

    for index, line in enumerate(lines):
        if any(keyword in line for keyword in section_keywords):
            for i in range(2,16):
                parts = lines[index+i].strip().split('\t')
                if len(parts) >= 3:
                    gene = parts[1]
                    site = parts[2]
                    if (gene,site) not in gene_site_pairs:
                        gene_site_pairs.append((gene, site))
                    else:
                        continue

    return gene_site_pairs

def extract_sequence_from_pdb(pdb_file):
    """
    Extract the amino acid sequence from a PDB file.

    Parameters:
    - pdb_file (str): Path to the input PDB file.

    Returns:
    - str: Amino acid sequence.
    """
    sequence = ""
    with open(pdb_file, "r") as pdb:
        for line in pdb:
            if line.startswith("ATOM") and "CA" in line:
                # Parse the line to extract the amino acid code
                fields = line.split()
                if len(fields) >= 4:
                    amino_acid = fields[3]
                    if amino_acid in three_to_one_letter:
                        sequence += three_to_one_letter[amino_acid]
    return sequence

def calculate_rmsd_bio(coords1, coords2):
    """
    Calculates the RMSD between two sets of coordinates using Biopython's Superimposer.
    """
    sup = Superimposer()
    sup.set_atoms(coords1, coords2)
    return sup.rms

def main():
    # Dictionnaire de protéines et de leurs sites
    protein_sites = {
        "glycogene_phosphorylase": "14",
    }

    with open(r"/mnt/c/Users/crisd/Desktop/rmsd_comparisons.txt", "w") as f:
        for protein, site in protein_sites.items():
            protein_folder = f"/mnt/c/Users/crisd/Desktop/Protein_comparaison/{protein}"

            if protein == "ERK1":
                real_unphospho_file = os.path.join(protein_folder, f"real-{protein}.pdb")
                real_phospho_file = os.path.join(protein_folder, f"real-phospho-{protein}-{site}.pdb")
                af_unphospho_file = os.path.join(protein_folder, f"{protein}-AF.pdb")
                af_d_phospho_file = os.path.join(protein_folder, f"ERK1-D183-D185.pdb")
                af_e_phospho_file = os.path.join(protein_folder, f"ERK1-E183-E185.pdb")
                contigs="A1-"+str(b)+"/A"+str(c)+"-"+str(d)
            else:
                real_unphospho_file = os.path.join(protein_folder, f"real-{protein}.pdb")
                real_phospho_file = os.path.join(protein_folder, f"real-phospho-{protein}-{site}.pdb")
                af_unphospho_file = os.path.join(protein_folder, f"{protein}-AF.pdb")
                af_d_phospho_file = os.path.join(protein_folder, f"{protein}-D{site}.pdb")
                af_e_phospho_file = os.path.join(protein_folder, f"{protein}-E{site}.pdb")
                sequence = extract_residues_from_PDB(af_d_phospho_file)
                d=len(sequence)
                print(d)
                print(sequence)
                site=int(site)
                b=site-1
                c=site+1
                contigs="A1-"+str(b)+"/A"+str(c)+"-"+str(d)
            print(contigs)
            f.write(f"Protein {protein}, phosphosite {site}:\n")
            
            # AF unphosphorylated versus AF phosphomimetic structure D
            coords1 = extract_sequence_from_pdb(af_unphospho_file)
            coords2 = extract_sequence_from_pdb(af_d_phospho_file)
            rmsd = calculate_rmsd_bio(coords1, coords2)
            f.write(f"\t AF unphosphorylated structure versus AF phosphomimetic structure D: {rmsd}\n")

            # AF unphosphorylated versus AF phosphomimetic structure E
            coords2 = extract_sequence_from_pdb(af_e_phospho_file)
            rmsd = calculate_rmsd_bio(coords1, coords2)
            f.write(f"\t AF unphosphorylated structure versus AF phosphomimetic structure E: {rmsd}\n")

            f.write("\n")

if __name__ == "__main__":
    main()
