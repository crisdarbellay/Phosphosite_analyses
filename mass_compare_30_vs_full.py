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

def sort_key_with_nan(value):
    if isinstance(value, float) and math.isnan(value):
        return (0, float('inf'))  # Place nan values at the end
    return (1, value)

def main():
    input_file = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/3_algorithm/results/filtered_results.txt"
    temp_file = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/3_algorithm/results/temp_results.txt"
    output_file = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_rmsd.txt"
    output_total = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_total.txt"
    output_file_sorted_d = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_rmsd_sorted_d.txt"
    output_file_sorted_e = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_rmsd_sorted_e.txt"
    output_txt = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_table.txt"
    output_csv = "/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_table.csv"

    with open(input_file, 'r') as f_in, open(temp_file, 'w') as f_temp, open(output_total, 'w') as tot:
        lines = f_in.readlines()
        header = lines[0].strip().split('\t')
        header.extend([
            "wt30_full_30", "wt30-D30_30", "wt30-E30_30",
            "wt30-D30_5", "wt30-E30_5"
        ])

        f_temp.write('\t'.join(header) + '\n')
        tot.write('\t'.join(header) + '\n')

        results = []

        for line in lines[1:]:
            if line.startswith("Kinase"):
                f_temp.write(line.split()[0])
                continue
            gene_data = line.strip().split('\t')
            gene_name = gene_data[1]
            site = gene_data[2]

            # Paths to PDB files for each type of data
            if os.path.exists( f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final/pdb/{gene_name}-{site}"):
                wild_type_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final/pdb/{gene_name}-{site}"
            else:
                wild_type_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final/pdb/{gene_name}"
            phosphomimetic_d_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final/pdb/{gene_name}-D{site}"
            phosphomimetic_e_folder = f"/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/2_algorithm/mass_test_final/pdb/{gene_name}-E{site}"
            wild_type_30_folder = f"/mnt/c/Users/crisd/Desktop/little-adapted/{gene_name}-{site}"
            phosphomimetic_d_30_folder = f"/mnt/c/Users/crisd/Desktop/little-adapted/{gene_name}-D{site}"
            phosphomimetic_e_30_folder = f"/mnt/c/Users/crisd/Desktop/little-adapted/{gene_name}-E{site}"

            if not os.path.exists(wild_type_folder):
                continue
            if not os.path.exists(phosphomimetic_d_folder):
                continue
            if not os.path.exists(phosphomimetic_e_folder):
                continue
            if not os.path.exists(wild_type_30_folder):
                continue
            if not os.path.exists(phosphomimetic_d_30_folder):
                continue
            if not os.path.exists(phosphomimetic_e_30_folder):
                continue

            # Fetch the appropriate PDB files for comparison
            def fetch_pdb(folder):
                for file in os.listdir(folder):
                    if "_relaxed_rank_001_"  in file and file.endswith(".pdb") and "sorted" not in file and "no_oxt" in file:
                        return os.path.join(folder, file)
                return None

            wild_type_pdb_full = fetch_pdb(wild_type_folder)
            phosphomimetic_d_pdb_full = fetch_pdb(phosphomimetic_d_folder)
            phosphomimetic_e_pdb_full = fetch_pdb(phosphomimetic_e_folder)

            wild_type_pdb_30 = fetch_pdb(wild_type_30_folder)
            phosphomimetic_d_pdb_30 = fetch_pdb(phosphomimetic_d_30_folder)
            phosphomimetic_e_pdb_30 = fetch_pdb(phosphomimetic_e_30_folder)

            if wild_type_pdb_full and phosphomimetic_d_pdb_full and phosphomimetic_e_pdb_full:
                sequence = extract_residues_from_PDB(wild_type_pdb_full)
                site = int(site)
                d = int(len(sequence))
                a = max(site - 5, 1)
                b = min(site + 5, d)
                e = max(site - 29, 1)
                f = min(site + 29, d)
                contigs_full = f"A1-{site - 1}/A{site + 1}-{d}"
                contigs_30 = f"A{e}-{site - 1}/A{site + 1}-{f}"
                contigs_5 = f"A{a}-{site - 1}/A{site + 1}-{b}"

                # RMSD comparisons between full pdb and different contigs
                try:
                    rmsd_d_full_contigs = round(CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                except Exception as e:
                    rmsd_d_full_contigs = 1000

                try:
                    rmsd_e_full_contigs = round(CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)
                except Exception as e:
                    rmsd_e_full_contigs = 1000

                try:
                    rmsd_d_30 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                except Exception as e:
                    rmsd_d_30 = 1000

                try:
                    rmsd_e_30 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)
                except Exception as e:
                    rmsd_e_30 = 1000

                try:
                    rmsd_d_5 = round(CompareTwoPDBs(contigs_5, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                except Exception as e:
                    rmsd_d_5 = 1000

                try:
                    rmsd_e_5 = round(CompareTwoPDBs(contigs_5, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)
                except Exception as e:
                    rmsd_e_5 = 1000

                # Additional RMSD comparisons with pdb_30
                try:
                    rmsd_wt30_full_contigs30 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_full, wild_type_pdb_30), 3)
                except Exception as e:
                    rmsd_wt30_full_contigs30 = 1000

                try:
                    rmsd_wt30_D30_contigs30 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_30, phosphomimetic_d_pdb_30), 3)
                except Exception as e:
                    rmsd_wt30_D30_contigs30 = 1000

                try:
                    rmsd_wt30_E30_contigs30 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_30, phosphomimetic_e_pdb_30), 3)
                except Exception as e:
                    rmsd_wt30_E30_contigs30 = 1000

                try:
                    rmsd_wt30_D30_contigs5 = round(CompareTwoPDBs(contigs_5, wild_type_pdb_30, phosphomimetic_d_pdb_30), 3)
                except Exception as e:
                    rmsd_wt30_D30_contigs5 = 1000

                try:
                    rmsd_wt30_E30_contigs5 = round(CompareTwoPDBs(contigs_5, wild_type_pdb_30, phosphomimetic_e_pdb_30), 3)
                except Exception as e:
                    rmsd_wt30_E30_contigs5 = 1000

                gene_data.extend([
                    rmsd_wt30_full_contigs30, rmsd_wt30_D30_contigs30, rmsd_wt30_E30_contigs30,
                    rmsd_wt30_D30_contigs5, rmsd_wt30_E30_contigs5
                ])
                f_temp.write('\t'.join(map(str, gene_data)) + '\n')
                results.append(gene_data)
                print(f"{gene_name} done")
            else:
                print(f"Missing PDB files for gene: {gene_name}")

    os.rename(temp_file, output_file)
    
    # Sorting results according to RMSD D and E (full)
    results_sorted_d = sorted(results, key=lambda x: sort_key_with_nan(x[-4]), reverse=True)
    write_sorted_results(output_file_sorted_d, header, results_sorted_d, use_pretty_table=True)

    results_sorted_e = sorted(results, key=lambda x: sort_key_with_nan(x[-3]), reverse=True)
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


if __name__ == "__main__":
    main()
