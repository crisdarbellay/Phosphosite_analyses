import os
from id_given_several_kinases import utility_id_given
import subprocess

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

def collabfold(sequence):
    tmp_folder = "collabfold_tmp"
    os.makedirs(tmp_folder, exist_ok=True)

    with open(os.path.join(tmp_folder, "temp_sequence.fa"), "w") as seq_file:
        seq_file.write(f">TempSequence\n{sequence}\n")

    subprocess.run(
        [
            "colabfold_batch",
            "--amber",
            "--rank",
            "plddt",
            "--use-gpu-relax",
            "--overwrite-existing-results",
            "--num-models",
            "1",
            "--num-relax",
            "1",
            os.path.join(tmp_folder, "temp_sequence.fa"),
            os.path.join(tmp_folder, "temp_results"),
        ]
    )

    pdb_file_path = os.path.join(tmp_folder, "temp_results", "temp_sequence_0_CF.pdb")

    return pdb_file_path
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
            if line.startswith("ATOM"):
                # Parse the line to extract the amino acid code
                fields = line.split()
                if len(fields) >= 4:
                    amino_acid = fields[3]
                    if amino_acid in three_to_one_letter:
                        sequence += three_to_one_letter[amino_acid]

    return sequence
# Function to perform a phosphomimetic
def phosphomimetic(pdb_file, site):
    """
    Perform phosphomimetic simulation on the given PDB file at the specified site.

    Parameters:
    - pdb_file (str): Path to the input PDB file.
    - site (int): Site of phosphorylation.

    Returns:
    - str: Path to the phosphomimetic PDB file.
    """
    sequence= extract_sequence_from_pdb(pdb_file)
    sequence[site]='D'
    result_pdb = collabfold(sequence)
    return result_pdb

def extract_gene_site_pairs(file_path):

    with open(file_path, 'r') as file:
        lines = file.readlines()
    gene_site_pairs = []

    pairs_extracted = 0
    section_keywords = ["Ranking Score", "Ranking AFConf", "Ranking NextAA"]

    for index, line in enumerate(lines):

        if any(keyword in line for keyword in section_keywords):
            for i in range(1,16):
                parts = lines[index+i].strip().split('\t')
                if len(parts) >= 3:
                    gene = parts[1]
                    site = parts[2]
                    if (gene,site) not in gene_site_pairs:
                        gene_site_pairs.append((gene, site))
                    else:
                        continue

    return gene_site_pairs

def save_cut_pdb(cut_pdb, outside_folder, protein, site):
    """
    Save the content of a PDB into a file in the specified outside folder.

    Parameters:
    - cut_pdb (str): Content of the PDB.
    - outside_folder (str): Path to the outside folder where the PDB will be saved.
    - protein (str): Name of the protein.
    - site (int): Site of phosphorylation.

    Returns:
    - None
    """
    # Create the outside folder if it doesn't exist
    os.makedirs(outside_folder, exist_ok=True)
    
    # Define the output file path
    output_file = os.path.join(outside_folder, f"{protein}_site{site}_cut.pdb")
    
    # Save the PDB content to the output file
    with open(output_file, "w") as output_pdb:
        output_pdb.write(cut_pdb)

# Function to find the PDB file for a protein
def find_pdb_from_protein(protein):
    """
    Find the corresponding PDB file for a given protein.

    Parameters:
    - protein (str): Name of the protein.

    Returns:
    - str: Path to the PDB file.
    """
    pdb_file_path = os.path.join("pdb_folder", f"{protein}.pdb")
    return pdb_file_path

# Function to print an error message to the output file
def print_outside_error(pdb_file, rmsd_un):
    """
    Print an error message to the output file.

    Parameters:
    - pdb_file (str): Path to the PDB file.
    - rmsd_un (float): RMSD value.

    Returns:
    - None
    """
    with open("result.txt", "a") as result_file:
        result_file.write(f"Error for file {pdb_file}, rmsd={rmsd_un}, too low\n")

# Function to cut the protein around the phosphorylation site
def cut_ball_around(pdb_file, site, size):
    """
    Cut the protein structure around the phosphorylation site.

    Parameters:
    - pdb_file (str): Path to the input PDB file.
    - site (int): Site of phosphorylation.
    - size (int): Size for cutting around the site.

    Returns:
    - str: Path to the cut PDB file.
    """
    sequence=extract_sequence_from_pdb(pdb_file)
    try:
        cut_sequence=sequence[site-size-1,site+size-1]
    except:
        try:
            cut_sequence=sequence[0,site+size-1]
        except:
            cut_sequence=sequence[site-size-1:]
    cut_pdb_file = collabfold(cut_sequence)
    return cut_pdb_file

# Main loop
x=4
with open("result.txt", "w") as result_file:
    for protein, site in top_15:
        original_pdb = find_pdb_from_protein(protein)
        a = 10  # Initialize a
        while True:
            temp_phospho = phosphomimetic(original_pdb, site)
            rmsd_un = compare_two_pdb(original_pdb, temp_phospho, contigs=((site - a), (site + a)))
            if a > 50:
                print_outside_error(original_pdb, rmsd_un)
                break
            if rmsd_un >= x:
                cut_pdb = cut_ball_around(original_pdb, site, a)
                cut_phospho = phosphomimetic(cut_pdb, site)
                rmsd_deux = compare_two_pdb(cut_pdb, cut_phospho)
                if rmsd_deux > x:
                    result_file.write(f"Positive result for protein {protein}, site {site}, rmsd_un={rmsd_un}, rmsd_deux={rmsd_deux}, a={a}\n")
                    save_cut_pdb(cut_pdb)
                else:
                    result_file.write(f"Negative result for protein {protein}, site {site}, rmsd_un={rmsd_un}, rmsd_deux={rmsd_deux}, a={a}\n")
                break
            a += 10
