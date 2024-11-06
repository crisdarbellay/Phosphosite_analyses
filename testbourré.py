import os
import tempfile
import subprocess
import csv
from Bio.PDB import PDBParser
from prettytable import PrettyTable
from main_tools2 import CompareTwoPDBs
from main_tools import extract_residues_from_PDB
from id_given_several_kinases import utility_id_given as util
import math

sec_struct_distance = 6  # Distance for secondary structure detection
count = 0
secondary_structures = [('H', sec_struct_distance, 'alpha', count), ('E', sec_struct_distance, 'beta', count), 
                        ('B', sec_struct_distance, 'iso_B', count), ('G', sec_struct_distance, 'alpha3', count), 
                        ('I', sec_struct_distance, 'alphaI', count), ('T', sec_struct_distance, 'hydturn', count)]

def sort_key_with_nan(value):
    if isinstance(value, float) and math.isnan(value):
        return (0, float('inf'))  # Place nan values at the end
    return (1, value)

def calculate_stability_pdb(pdb_file_path, site_data, secondary_structures):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file_path)
    model = structure[0]
    stability_scores = []

    temp_dssp_file_path = tempfile.mkstemp()[1]
    run_dssp(pdb_file_path, temp_dssp_file_path)
    
    dssp_info = parse_dssp(temp_dssp_file_path)
    stability_score, chain_distance, average_confidence_score = calculate_scores(model, site_data, secondary_structures, dssp_info, pdb_file_path)
    
    stability_scores.append({
        'residue_number': site_data,
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

def calculate_scores(model, site_data, secondary_structures, dssp_info, pdb_file_path):
    # Vérifie si site_data est une liste
    if isinstance(site_data, list):
        site_residue_number = int(site_data[0]["sub_site"])  # Utilise le premier élément
    else:
        site_residue_number = int(site_data)
    stability_score = 0
    chain_distance = 100000
    confidence_scores_sum = 0
    confidence_scores_count = 0
    count_dict = {count[2]: 0 for count in secondary_structures}
    
    residues = model.get_residues()
    
    for residu1 in residues:
        if residu1.get_id()[1] == site_residue_number:
            for residu2 in model.get_residues():
                distance = util.calculate_distance(residu1['CA'].get_coord(), residu2['CA'].get_coord())
                res_dssp = dssp_info.get(residu2.get_id()[1], '')
                
                for count in secondary_structures:
                    if count[0] == res_dssp and distance <= 6:
                        chain_distance_temp = abs(residu1.get_id()[1] - residu2.get_id()[1])
                        if chain_distance_temp < chain_distance:
                            chain_distance = chain_distance_temp
                        count_dict[count[2]] += 1

                if distance <= 10.0:
                    confidence_score = util.calculate_confidence_scores(util.begin_by_ATOM(pdb_file_path), residu2.get_id())
                    if confidence_score is not None:
                        confidence_scores_sum += confidence_score
                        confidence_scores_count += 1

    stability_score = sum(count_dict.values())
    average_confidence_score = round(confidence_scores_sum / confidence_scores_count, 1) if confidence_scores_count > 0 else None
    
    return stability_score, chain_distance, average_confidence_score

def main():
    input_file = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/selection_no_doublon.txt"
    temp_file = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/phosphosites/results/temp_results.txt"
    output_file = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/phosphosites/results/results_rmsd_score.txt"
    output_total = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/phosphosites/results/results_total__score.txt"
    output_file_sorted_d = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/phosphosites/results/results_rmsd_sorted_d_score.txt"
    output_file_sorted_e = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/phosphosites/results/results_rmsd_sorted_e_score.txt"
    output_txt = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/phosphosites/results/results_table.txt"
    output_csv = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/phosphosites/results/results_table.csv"

    with open(input_file, 'r') as f_in, open(temp_file, 'w') as f_temp, open(output_total, 'w') as tot:
        lines = f_in.readlines()
        header = lines[0].strip().split('\t')
        header.extend([
        "wt-D_full", "wt-E_full", "wt-D_30", "wt-E_30", "wt-D_5", "wt-E_5",
        "score_diff_D", "score_diff_E"  
        ])
        f_temp.write('\t'.join(header) + '\n')
        tot.write('\t'.join(header) + '\n')

        results = []

        for line in lines[1:]:
            kinase=line.split()[0]
            if kinase != 'Src' and kinase != 'PKACA':
                continue
            if line.startswith("Kinase"):
                f_temp.write(line.split()[0])
                continue
            gene_data = line.strip().split('\t')
            gene_name = gene_data[1]
            site = gene_data[2]
            score = gene_data[3]

            # Paths to PDB files for each type of data
            wild_type_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final_2/pdb/{gene_name}"
            phosphomimetic_d_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final_2/pdb/{gene_name}-D{site}"
            phosphomimetic_e_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final_2/pdb/{gene_name}-E{site}"

            if not os.path.exists(wild_type_folder):
                continue
            if not os.path.exists(phosphomimetic_d_folder):
                continue
            if not os.path.exists(phosphomimetic_e_folder):
                continue

            # Fetch the appropriate PDB files for comparison
            def fetch_pdb(folder):
                for file in os.listdir(folder):
                    if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                        return os.path.join(folder, file)
                return None

            wild_type_pdb_full = fetch_pdb(wild_type_folder)
            phosphomimetic_d_pdb_full = fetch_pdb(phosphomimetic_d_folder)
            phosphomimetic_e_pdb_full = fetch_pdb(phosphomimetic_e_folder)

            if wild_type_pdb_full and phosphomimetic_d_pdb_full and phosphomimetic_e_pdb_full:
                sequence = extract_residues_from_PDB(wild_type_pdb_full)
                site = int(site)
                d = int(len(sequence))
                a = max(site - 5, 1)
                b = min(site + 5, d)
                e = max(site - 30, 1)
                f = min(site + 30, d)
                contigs_full = f"A1-{site - 1}/A{site +         1}-{d}"
                contigs_30 = f"A{e}-{site - 1}/A{site + 1}-{f}"
                contigs_5 = f"A{a}-{site - 1}/A{site + 1}-{b}"

                # RMSD comparisons between full pdb and different contigs
                rmsd_d_full_contigs = round(CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                rmsd_e_full_contigs = round(CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)
                rmsd_d_30 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                rmsd_e_30 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)
                rmsd_d_5 = round(CompareTwoPDBs(contigs_5, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                rmsd_e_5 = round(CompareTwoPDBs(contigs_5, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)

                # New score calculations using calculate_stability_pdb function
                score_D = calculate_stability_pdb(phosphomimetic_d_pdb_full, site, secondary_structures)
                score_E = calculate_stability_pdb(phosphomimetic_e_pdb_full, site, secondary_structures)
                score_diff_D = float(score) - float(score_D[0]['stability_score'])
                score_diff_E = float(score) - float(score_E[0]['stability_score'])

                gene_data.extend([
                    rmsd_d_full_contigs, rmsd_e_full_contigs, rmsd_d_30, rmsd_e_30, rmsd_d_5, rmsd_e_5, score_diff_D, score_diff_E
                ])
                f_temp.write('\t'.join(map(str, gene_data)) + '\n')
                results.append(gene_data)
                print(f"{gene_name} done")
            else:
                print(f"Missing PDB files for gene: {gene_name}")

    os.rename(temp_file, output_file)
    
    # Sorting results according to RMSD D and E (full)
    results_sorted_d = sorted(results, key=lambda x: sort_key_with_nan(x[-5]), reverse=True)
    write_sorted_results(output_file_sorted_d, header, results_sorted_d, use_pretty_table=True)

    results_sorted_e = sorted(results, key=lambda x: sort_key_with_nan(x[-4]), reverse=True)
    write_sorted_results(output_file_sorted_e, header, results_sorted_e, use_pretty_table=True)

def write_sorted_results(file_path, header, results_sorted, use_pretty_table=False):
    if use_pretty_table:
        table = PrettyTable()
        table.field_names = header
        for result in results_sorted:
            table.add_row(result)
        table.border = False
        table.header = True
        with open(file_path, 'w') as f_out_sorted:
            f_out_sorted.write(table.get_string())
    else:
        with open(file_path, 'w') as f_out_sorted:
            f_out_sorted.write('\t'.join(header) + '\n')
            for result in results_sorted:
                f_out_sorted.write('\t'.join(map(str, result)) + '\n')

