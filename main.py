import tempfile
import math
import subprocess
import os
from Bio.PDB import PDBParser



def begin_by_ATOM(pdb_path):
    """
    Extract the lines starting with 'ATOM' from a PDB file.
    """
    atom_lines = []
    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM'):
                atom_lines.append(line)
    return atom_lines

def calculate_distance(coord1, coord2):
    """
    Calculate the Euclidean distance between two 3D coordinates.
    """
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    return distance

def extract_confidence_score(line):
    """
    Extract the confidence score from a line in the PDB file.
    """
    confidence_score = float(line[61:66])
    return confidence_score

def calculate_confidence_scores(atom_lines, residue_key):
    """
    Calculate the confidence scores for a specific residue.
    """
    confidence_scores = {}
    for line in atom_lines:
        line_residue_number = int(line.split()[5])
        if line_residue_number == residue_key[1]:
            confidence_score = extract_confidence_score(line)
            confidence_scores[residue_key] = confidence_score
            break
    else:
        confidence_scores[residue_key] = 0.0
    return confidence_scores

def calculate_stability(pdb_path, location_list):
    """
    Calculate the stability scores for residues at specific locations in a PDB file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]  # Assuming there is only one model in the structure

    stability_scores = []

    # Generate DSSP file
    temp_dssp_file_path = tempfile.mkstemp()[1]

    if os.name == 'nt':
        # Windows environment
        wsl_dssp_command = f"wsl mkdssp {pdb_path} {temp_dssp_file_path}"
        subprocess.call(wsl_dssp_command, shell=True)
    else:
        # Non-Windows environment
        dssp_executable = "mkdssp"
        subprocess.call([dssp_executable, pdb_path, temp_dssp_file_path])

    for residue_number in location_list:
        residue_key = (' ', residue_number, ' ')

        helix_count = 0
        beta_count = 0
        stability_score = 0

        # Calculate distance cutoff for helix and beta sheet
        alpha_helix_distance = 10.0
        beta_fold_distance = 10.0

        # Get residue from model
        residues = model.get_residues()
        for residu1 in residues:
            if residu1.get_id()[1] == residue_number:
                for residu2 in model.get_residues():
                    distance = calculate_distance(residu1['CA'].get_coord(), residu2['CA'].get_coord())
                    if distance <= alpha_helix_distance:
                        res_dssp = ''
                        with open(temp_dssp_file_path, 'r') as dssp_file:
                            for dssp_line in dssp_file:
                                if dssp_line[5:10].strip() == str(residu2.get_id()[1]):
                                    res_dssp = dssp_line.split()[4]
                                    break
                        if res_dssp == 'H':
                            helix_count += 1
                    if distance <= beta_fold_distance:
                        res_dssp = ''
                        with open(temp_dssp_file_path, 'r') as dssp_file:
                            for dssp_line in dssp_file:
                                if dssp_line[5:10].strip() == str(residu2.get_id()[1]):
                                    res_dssp = dssp_line.split()[4]
                                    break
                        if res_dssp == 'E':
                            beta_count += 1

        stability_score = helix_count * 10 + beta_count * 5

        stability_scores.append({
            'residue_number': residue_number,
            'stability_score': stability_score,
            'helix_count': helix_count,
            'beta_count': beta_count,
            'confidence_score': calculate_confidence_scores(begin_by_ATOM(pdb_path), residue_key)[residue_key]
        })

    # Remove temporary DSSP file
    os.remove(temp_dssp_file_path)

    return stability_scores

# Example usage
pdb_path = r"/mnt/c/Users/crisd/Desktop/Projet/NTH3.pdb"

location_list = [27, 90,134, 150, 300,190,246]  # Add more residue positions as needed

result = calculate_stability(pdb_path, location_list)

for stats in result:
    residue_number = stats['residue_number']
    print(f"Residue {residue_number}:")
    print(f"Stability Score: {stats['stability_score']}")
    print(f"Helix Count: {stats['helix_count']}")
    print(f"Beta Count: {stats['beta_count']}")
    print(f"Confidence Score: {stats['confidence_score']}")
    print("\n")