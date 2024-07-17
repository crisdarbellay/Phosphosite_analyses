"""
import os

def add_header_to_files(directory):
    # Parcourir tous les fichiers et sous-répertoires du répertoire
    for root, dirs, files in os.walk(directory):
        root_name = os.path.basename(root)
        if "AATF" in root_name:
            continue
        header = f"HEADER    DNA BINDING PROTEIN/DNA                 15-MAY-24   {root_name}"
        for filename in files:
            if ".pdb" not in filename:
                continue
            file_path = os.path.join(root, filename)
            
            # Vérifier si le chemin est un fichier
            if os.path.isfile(file_path):
                # Lire le contenu original du fichier
                with open(file_path, 'r') as file:
                    original_content = file.read()
                
                # Ajouter l'en-tête au début du contenu original
                new_content = header + "\n" + original_content
                
                # Écrire le nouveau contenu dans le fichier
                with open(file_path, 'w') as file:
                    file.write(new_content)
                print(f"Header added to {file_path}")

if __name__ == "__main__":
    directory_path = "/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-30"
    add_header_to_files(directory_path)
"""

import os
import tempfile
import subprocess
from Bio.PDB import PDBParser
from prettytable import PrettyTable
from main_tools2 import CompareTwoPDBs
from main_tools import extract_residues_from_PDB
from id_given_several_kinases import utility_id_given as util

# Constants
THREE_TO_ONE_LETTER = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

SEC_STRUCT_DISTANCE = 8
SECONDARY_STRUCTURES = [
    ('H', SEC_STRUCT_DISTANCE, 'alpha', 0),
    ('E', SEC_STRUCT_DISTANCE, 'beta', 0),
    ('B', SEC_STRUCT_DISTANCE, 'iso_B', 0),
    ('G', SEC_STRUCT_DISTANCE, 'alpha3', 0),
    ('I', SEC_STRUCT_DISTANCE, 'alphaI', 0),
    ('T', SEC_STRUCT_DISTANCE, 'hydturn', 0)
]

# Functions

def calculate_stability_pdb(pdb_file_path, site_data, secondary_structures):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file_path)
    model = structure[0]
    stability_scores = []

    temp_dssp_file_path = tempfile.mkstemp()[1]
    run_dssp(pdb_file_path, temp_dssp_file_path)
    
    dssp_info = parse_dssp(temp_dssp_file_path)
    stability_score, chain_distance, average_confidence_score = calculate_scores(model, site_data, secondary_structures, dssp_info)
    
    stability_scores.append({
        'residue_number': site_data[1],
        'stability_score': stability_score,
        'average_confidence_score': average_confidence_score,
        'chain_dist_secondary_struct': chain_distance,
        **{count[2]: 0 for count in secondary_structures}
    })
    
    os.remove(temp_dssp_file_path)
    return stability_scores

def run_dssp(pdb_file_path, output_path):
    if os.name == 'nt':
        wsl_dssp_command = f"wsl mkdssp {pdb_file_path} {output_path}"
        subprocess.call(wsl_dssp_command, shell=True)
    else:
        subprocess.call(["mkdssp", pdb_file_path, output_path, "--output-format=dssp"])

def parse_dssp(dssp_file_path):
    dssp_info = {}
    with open(dssp_file_path, 'r') as dssp_file:
        dssp_lines = dssp_file.readlines()
        for line in dssp_lines:
            parts = line.split()
            if len(parts) > 4 and parts[1].isdigit():
                residue_number = int(parts[1])
                res_dssp = parts[4]
                dssp_info[residue_number] = res_dssp
    return dssp_info

def calculate_scores(model, site_data, secondary_structures, dssp_info):
    stability_score = 0
    chain_distance = 100000
    confidence_scores_sum = 0
    confidence_scores_count = 0
    count_dict = {count[2]: 0 for count in secondary_structures}
    
    residues = model.get_residues()
    site_residue_number = int(site_data[1])
    
    for residu1 in residues:
        if residu1.get_id()[1] == site_residue_number:
            for residu2 in model.get_residues():
                distance = util.calculate_distance(residu1['CA'].get_coord(), residu2['CA'].get_coord())
                res_dssp = dssp_info.get(residu2.get_id()[1], '')
                
                for count in secondary_structures:
                    if count[0] == res_dssp and distance <= count[1]:
                        chain_distance_temp = abs(residu1.get_id()[1] - residu2.get_id()[1])
                        if chain_distance_temp < chain_distance:
                            chain_distance = chain_distance_temp
                        count_dict[count[2]] += 1

                if distance <= 10.0:
                    confidence_score = util.calculate_confidence_scores(util.begin_by_ATOM(pdb_file_path), residu2.get_id()).get(residu2.get_id(), None)
                    if confidence_score is not None:
                        confidence_scores_sum += confidence_score
                        confidence_scores_count += 1

    stability_score = sum(count_dict.values())
    average_confidence_score = round(confidence_scores_sum / confidence_scores_count, 1) if confidence_scores_count > 0 else None
    
    return stability_score, chain_distance, average_confidence_score

def process_kinase_data(line, gene_data, secondary_structures):
    gene_name = gene_data[0]
    site = gene_data[1]
    folders = prepare_folders(gene_name, site)
    
    pdb_files = find_pdb_files(folders)
    if not all(pdb_files.values()):
        print(f"Missing PDB files for gene: {gene_name}")
        return None
    
    rmsd_results = calculate_rmsd(gene_data, pdb_files)
    stability_results = calculate_stability_scores(pdb_files, gene_data, secondary_structures)
    
    gene_data.extend(rmsd_results + stability_results)
    return gene_data

def prepare_folders(gene_name, site):
    base_path = "/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb"
    adapted_path_30 = "/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-30-contigs-adapted"
    adapted_path_50 = "/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-50-contigs-adapted"
    
    return {
        "wt": os.path.join(base_path, gene_name),
        "d": os.path.join(base_path, f"{gene_name}-D{site}"),
        "e": os.path.join(base_path, f"{gene_name}-E{site}"),
        "wt_30": os.path.join(adapted_path_30, f"{gene_name}-{site}"),
        "d_30": os.path.join(adapted_path_30, f"{gene_name}-D{site}"),
        "e_30": os.path.join(adapted_path_30, f"{gene_name}-E{site}"),
        "wt_50": os.path.join(adapted_path_50, f"{gene_name}-{site}"),
        "d_50": os.path.join(adapted_path_50, f"{gene_name}-D{site}"),
        "e_50": os.path.join(adapted_path_50, f"{gene_name}-E{site}")
    }

def find_pdb_files(folders):
    pdb_files = {}
    for key, folder in folders.items():
        if not os.path.exists(folder):
            pdb_files[key] = None
        else:
            for file in os.listdir(folder):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    pdb_files[key] = os.path.join(folder, file)
                    break
    return pdb_files

def calculate_rmsd(gene_data, pdb_files):
    sequence = extract_residues_from_PDB(pdb_files["wt"])
    site = int(gene_data[1])
    d = len(sequence)
    a, b = max(site - 30, 1), min(site + 29, d)
    e, f = max(site - 49, 1), min(site + 49, d)
    
    contigs_full = f"A1-{site - 1}/A{site + 1}-{d}"
    contigs_30 = f"A{a}-{site - 1}/A{site + 1}-{b}"
    contigs_50 = f"A{e}-{site - 1}/A{site + 1}-{f}"
    
    rmsd_results = []
    try:
        rmsd_results.extend([
            round(CompareTwoPDBs(contigs_full, pdb_files["wt"], pdb_files["d"]), 3),
            round(CompareTwoPDBs(contigs_full, pdb_files["wt"], pdb_files["e"]), 3),
            round(CompareTwoPDBs(contigs_30, pdb_files["wt"], pdb_files["d"]), 3),
            round(CompareTwoPDBs(contigs_30, pdb_files["wt"], pdb_files["e"]), 3)
        ])
    except Exception as e:
        print(f"Error calculating RMSD: {e}")
        rmsd_results.extend([0, 0, 0, 0])
    
    return rmsd_results

def calculate_stability_scores(pdb_files, gene_data, secondary_structures):
    stability_wt = calculate_stability_pdb(pdb_files["wt"], gene_data, secondary_structures)
    stability_d = calculate_stability_pdb(pdb_files["d"], gene_data, secondary_structures)
    stability_e = calculate_stability_pdb(pdb_files["e"], gene_data, secondary_structures)
    
    return [stability_wt, stability_d, stability_e, 0, 0]  # Replace 0, 0 with actual diff scores if needed

def main():
    input_file = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results.txt"
    temp_file = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/temp_results.txt"
    output_files = {
        "final": "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd.txt",
        "total": "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_total.txt",
        "sorted_d": "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd_sorted_d.txt",
        "sorted_e": "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd_sorted_e.txt",
        "sorted_d_diff": "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd_sorted_diff_d.txt",
        "sorted_e_diff": "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd_sorted_diff_e.txt"
    }
    
    with open(input_file, 'r') as f_in, open(temp_file, 'w') as f_temp, open(output_files["total"], 'w') as tot:
        lines = f_in.readlines()
        header = lines[0].strip().split('\t')
        header.extend([
            "scoreD", "scoreE", "rmsd-D", "rmsd-E", "D-30AA", "E-30AA", 
            "wt-d_f", "wt-e_f", "wt_f-wt_100", "wt_f-wt_60",
            "d_f-d_100", "d_f-d_60", "e_f-e_100", "e_f-e_60",
            "wt-d_100", "wt-e_100", "wt-d_60", "wt-e_60"
        ])
        
        f_temp.write('\t'.join(header) + '\n')
        tot.write('\t'.join(header) + '\n')
        
        results = []
        for line in lines[1:20]:
            if line.startswith("Kinase"):
                f_temp.write(line.split()[0])
                continue
            gene_data = line.strip().split('\t')
            processed_data = process_kinase_data(line, gene_data, SECONDARY_STRUCTURES)
            if processed_data:
                f_temp.write('\t'.join(map(str, processed_data)) + '\n')
                results.append(processed_data)
                print(f"{processed_data[0]} done")
        
        os.rename(temp_file, output_files["final"])
        sort_and_write_results(results, header, output_files)

def sort_and_write_results(results, header, output_files):
    results_sorted_d = sorted(results, key=lambda x: x[-22], reverse=True)
    write_sorted_results(output_files["sorted_d"], header, results_sorted_d)
    
    results_sorted_e = sorted(results, key=lambda x: x[-21], reverse=True)
    write_sorted_results(output_files["sorted_e"], header, results_sorted_e)
    
    results_sorted_d_diff = sorted(results, key=lambda x: x[-2], reverse=True)
    write_sorted_results(output_files["sorted_d_diff"], header, results_sorted_d_diff)
    
    results_sorted_e_diff = sorted(results, key=lambda x: x[-1], reverse=True)
    write_sorted_results(output_files["sorted_e_diff"], header, results_sorted_e_diff)

def write_sorted_results(file_path, header, results_sorted):
    table = PrettyTable()
    table.field_names = header
    for result in results_sorted:
        table.add_row(result)
    
    with open(file_path, 'w') as f_out_sorted:
        f_out_sorted.write(table.get_string())
        f_out_sorted.write('\n')

if __name__ == "__main__":
    main()
