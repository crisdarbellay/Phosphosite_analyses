import tempfile
import math
import subprocess
import os

from Bio.PDB import PDBParser
from Bio.PDB import PDBParser

def count_alpha_helices(pdb_path):
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_path)

    alpha_helix_count = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_id()[0].strip() == 'H_P':
                    alpha_helix_count += 1

    return alpha_helix_count


def begin_by_ATOM(pdb_path):
    atom_lines = []
    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM'):
                atom_lines.append(line)
    return atom_lines

def calculate_distance(coord1, coord2):
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    return distance

def extract_confidence_score(line):
    confidence_score = float(line[61:66])
    return confidence_score

def calculate_confidence_scores(atom_lines, residue_key):
    confidence_scores = {}
    for line in atom_lines:
        line_residue_number = int(line.split()[5])
        if line_residue_number == residue_key:
            confidence_score = extract_confidence_score(line)
            confidence_scores[residue_key] = confidence_score
            break
    else:
        confidence_scores[residue_key] = 0.0
    return confidence_scores

def calculate_stability(pdb_path, location_list):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]  # Assuming there is only one model in the structure

    stability_scores = {}

    print(f"Number of residues in the model: {len(list(model.get_residues()))}")
    print(f"Structure ID: {structure.id}")
    print(f"Number of models: {len(structure)}")
    print(f"Number of chains: {len(list(structure.get_chains()))}")
    print("\n")

    # Generate DSSP file
    temp_dssp_file = tempfile.NamedTemporaryFile(delete=False)
    temp_dssp_file_path = temp_dssp_file.name
    temp_dssp_file.close()

    if os.name == 'nt':
        # Windows environment
        wsl_dssp_command = f"wsl mkdssp {pdb_path} {temp_dssp_file_path}"
        subprocess.call(wsl_dssp_command, shell=True)
    else:
        # Non-Windows environment
        dssp_executable = "mkdssp"
        subprocess.call([dssp_executable, pdb_path, temp_dssp_file_path])

    for residue_number in location_list:
        residue_key = residue_number

        helix_count = 0
        beta_count = 0

        stability_score = 0

        # Calculate distance cutoff for helix and beta sheet
        distance_cutoff_helix = 1000.0
        distance_cutoff_beta = 1000.0

        # Get residue from model
        residue = None
        for chain in model:
            for res in chain:
                if res.get_id()[1] == residue_number:
                    residue = res
                    break
            if residue:
                break

        if residue is None:
            continue

        for target_residue_number in location_list:
            if residue_number != target_residue_number:
                target_residue = None
                for chain in model:
                    for res in chain:
                        if res.get_id()[1] == target_residue_number:
                            target_residue = res
                            break
                    if target_residue:
                        break

                if target_residue is None:
                    continue

                if residue.has_id("CA") and target_residue.has_id("CA"):
                    distance = calculate_distance(residue["CA"].coord, target_residue["CA"].coord)
                    if distance <= distance_cutoff_helix:
                        helix_count += 1
                        stability_score += 10

                    # Read DSSP attributes from the generated DSSP file
                    with open(temp_dssp_file_path, 'r') as dssp_file:
                        for dssp_line in dssp_file:
                            if dssp_line.startswith(f"  {residue_number} "):
                                dssp_attrs = dssp_line.split()
                                residue_dssp = dssp_attrs[2]
                                target_residue_dssp = dssp_attrs[2 + (target_residue_number - residue_number)]
                                if residue_dssp == "E" and target_residue_dssp == "E":
                                    if distance <= distance_cutoff_beta:
                                        beta_count += 1
                                        stability_score += 5

        # Read confidence score from PDB
        atom_lines = begin_by_ATOM(pdb_path)
        confidence_scores = calculate_confidence_scores(atom_lines, residue_key)
        confidence_score = confidence_scores[residue_key]

        stability_scores[residue_key] = {
            'stability_score': stability_score,
            'helix_count': helix_count,
            'beta_count': beta_count,
            'confidence_score': confidence_score
        }

    # Remove the temporary DSSP file
    os.remove(temp_dssp_file_path)

    return stability_scores


# Example usage
pdb_path = r"C:\Users\crisd\Desktop\Projet\NTH3.pdb"
location_list = [27, 90, 150, 300, 390]  # Add more residue positions as needed

result = calculate_stability(pdb_path, location_list)

for residue_number, stats in result.items():
    print(f"Residue {residue_number}:")
    print(f"Stability Score: {stats['stability_score']}")
    print(f"Helix Count: {stats['helix_count']}")
    print(f"Beta Count: {stats['beta_count']}")
    print(f"Confidence Score: {stats['confidence_score']}")
    print("\n")
