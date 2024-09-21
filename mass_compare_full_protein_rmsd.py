import os
import tempfile
import subprocess
import csv
from Bio.PDB import PDBParser
from prettytable import PrettyTable
from main_tools2 import CompareTwoPDBs
from main_tools import extract_residues_from_PDB
from id_given_several_kinases import utility_id_given as util
from no_id_given import utility_no_id_given as no

import math
sec_struct_distance = 6     #distance of detection for the secondary structures around your sites
count = 0
secondary_structures = [('H', sec_struct_distance, 'alpha', count), ('E', sec_struct_distance, 'beta', count), ('B', sec_struct_distance, 'iso_B', count), ('G', sec_struct_distance, 'alpha3', count), ('I', sec_struct_distance, 'alphaI', count), ('T', sec_struct_distance, 'hydturn', count)]

def sort_key_with_nan(value):
    if isinstance(value, float) and math.isnan(value):
        return (0, float('inf'))  # Place nan values at the end
    return (1, value)

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
            "wt-D_full", "wt-E_full", "wt-D_30", "wt-E_30", "wt-D_5", "wt-E_5"
        ])

        f_temp.write('\t'.join(header) + '\n')
        tot.write('\t'.join(header) + '\n')

        results = []

        for line in lines[1:20]:
            if line.startswith("Kinase"):
                f_temp.write(line.split()[0])
                continue
            gene_data = line.strip().split('\t')
            gene_name = gene_data[1]
            site = gene_data[2]
            score = gene_data[3]

            # Paths to PDB files for each type of data
            wild_type_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final/pdb/{gene_name}"
            phosphomimetic_d_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final/pdb/{gene_name}-D{site}"
            phosphomimetic_e_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final/pdb/{gene_name}-E{site}"

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
                contigs_full = f"A1-{site - 1}/A{site + 1}-{d}"
                contigs_30 = f"A{e}-{site - 1}/A{site + 1}-{f}"
                contigs_5 = f"A{a}-{site - 1}/A{site + 1}-{b}"

                # RMSD comparisons between full pdb and different contigs
                rmsd_d_full_contigs = round(CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                rmsd_e_full_contigs = round(CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)
                rmsd_d_30 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                rmsd_e_30 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)
                rmsd_d_5 = round(CompareTwoPDBs(contigs_5, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                rmsd_e_5 = round(CompareTwoPDBs(contigs_5, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)
                score_D = no.calculate_stability_cif(phosphomimetic_d_pdb_full, site,secondary_structures)
                score_E = no.calculate_stability_cif(phosphomimetic_e_pdb_full, site,secondary_structures)
                score_diff_D = score-score_D
                score_diff_E = score-score_D

                gene_data.extend([
                    rmsd_d_full_contigs, rmsd_e_full_contigs, rmsd_d_30, rmsd_e_30, rmsd_d_5, rmsd_e_5,score_diff_D,score_diff_E
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


main()
