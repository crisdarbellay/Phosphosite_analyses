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
sec_struct_distance = 8  # distance of detection for the secondary structures around your sites

count = 0
secondary_structures = [('H', sec_struct_distance, 'alpha', count), ('E', sec_struct_distance, 'beta', count),
                        ('B', sec_struct_distance, 'iso_B', count), ('G', sec_struct_distance, 'alpha3', count),
                        ('I', sec_struct_distance, 'alphaI', count), ('T', sec_struct_distance, 'hydturn', count)]


def calculate_stability_pdb(pdb_file_path, site_data, secondary_structures):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file_path)
    model = structure[0]  # Assuming there is only one model in the structure
    results = []

    stability_scores = []

    # Generate DSSP file
    temp_dssp_file_path = tempfile.mkstemp()[1]

    if os.name == 'nt':
        # Windows environment
        wsl_dssp_command = f"wsl mkdssp {cif_file_path} {temp_dssp_file_path}"
        subprocess.call(wsl_dssp_command, shell=True)
    else:
        # Non-Windows environment
        dssp_executable = "mkdssp"
        subprocess.call([dssp_executable, pdb_file_path, temp_dssp_file_path, "--output-format=dssp"])

    # Store DSSP information in a dictionary
    dssp_info = {}
    with open(temp_dssp_file_path, 'r') as dssp_file:
        dssp_lines = dssp_file.readlines()
        for dssp_line in dssp_lines:
            line = dssp_line.split()
            if line[1] == line[0] and line[2] == 'A':
                residue_number = int(line[1])
                res_dssp = line[4]
                dssp_info[residue_number] = res_dssp
    stability_score = 0
    chain_distance = 100000
    chain_distance_temp = 10000
    confidence_scores_sum = 0
    confidence_scores_count = 0
    count_dict = {count[2]: 0 for count in secondary_structures}
    residues = model.get_residues()
    if site_data[1] in [entry['residue_number'] for entry in stability_scores]:
        pass
    for residu1 in residues:
        resname1 = three_to_one_letter.get(residu1.get_resname(), '')
        if residu1.get_id()[1] == int(site_data[1]):
            for residu2 in model.get_residues():
                distance = util.calculate_distance(residu1['CA'].get_coord(), residu2['CA'].get_coord())
                res_dssp = dssp_info.get(residu2.get_id()[1], '')

                for count in secondary_structures:
                    if count[0] == res_dssp:
                        chain_distance_temp = abs(residu1.get_id()[1] - residu2.get_id()[1])
                        if chain_distance_temp < chain_distance:
                            chain_distance = chain_distance_temp
                    if distance <= count[1] and res_dssp == count[0]:
                        count_dict[count[2]] += 1
                if distance <= 10.0:
                    try:
                        confidence_score = util.calculate_confidence_scores(util.begin_by_ATOM(pdb_file_path),
                                                                            residu2.get_id())[residu2.get_id()]
                    except:
                        confidence_score = None
                    if confidence_score is not None:
                        confidence_scores_sum += confidence_score
                        confidence_scores_count += 1
    stability_score = sum(count_dict.values())
    average_confidence_score = round(confidence_scores_sum / confidence_scores_count,
                                     1) if confidence_scores_count > 0 else None

    stability_scores.append({
        'residue_number': site_data[1],
        'stability_score': stability_score,
        'average_confidence_score': average_confidence_score,
        'chain_dist_secondary_struct': chain_distance,
        **count_dict
    })
    results.append({
        'gene_name': site_data[0],
        'stability_scores': stability_scores,
    })
    # Remove temporary DSSP file
    os.remove(temp_dssp_file_path)

    return stability_scores


def main():
    input_file = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results.txt"
    temp_file = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/temp_results.txt"
    output_file = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd.txt"
    output_total = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_total.txt"
    output_file_sorted_d = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd_sorted_d.txt"
    output_file_sorted_e = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd_sorted_e.txt"
    output_file_sorted_d_diff = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd_sorted_wt-60.txt"
    output_file_sorted_e_diff = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_rmsd_sorted_wt-100.txt"
    output_txt = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_table.txt"
    output_csv = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final/results_table.csv"

    with open(input_file, 'r') as f_in, open(temp_file, 'w') as f_temp, open(output_total, 'w') as tot:
        lines = f_in.readlines()
        header = lines[0].strip().split('\t')
        header.extend([
            "wt-D_full", "wt-E_full", "wt_full-wt_100", "wt_full-wt_60", "D_full-D_100", "D_full-D_60", "E_full-E_100", "E_full-E_60", "wt_100-D_100",
             "wt_100-E_100", "wt_60-D_60", "wt_60-E_60"
        ])

                
        f_temp.write('\t'.join(header) + '\n')
        tot.write('\t'.join(header) + '\n')

        results = []

        for line in lines[1:]:
            if line.startswith("Kinase"):
                f_temp.write(line.split()[0])
                continue
            gene_data = line.strip().split('\t')
            gene_name = gene_data[0]
            site = gene_data[1]

            # Paths to PDB files for each type of data
            wild_type_folder = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb/{gene_name}"
            phosphomimetic_d_folder = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb/{gene_name}-D{site}"
            phosphomimetic_e_folder = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb/{gene_name}-E{site}"
            wild_type_folder_30 = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-30-contigs-adapted/{gene_name}-{site}"
            phosphomimetic_d_folder_30 = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-30-contigs-adapted/{gene_name}-D{site}"
            phosphomimetic_e_folder_30 = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-30-contigs-adapted/{gene_name}-E{site}"
            wild_type_folder_50 = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-50-contigs-adapted/{gene_name}-{site}"
            phosphomimetic_d_folder_50 = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-50-contigs-adapted/{gene_name}-D{site}"
            phosphomimetic_e_folder_50 = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-50-contigs-adapted/{gene_name}-E{site}"

            if not os.path.exists(wild_type_folder):
                continue
            if not os.path.exists(phosphomimetic_d_folder):
                continue
            if not os.path.exists(phosphomimetic_e_folder):
                continue
            if not os.path.exists(wild_type_folder_30):
                continue
            if not os.path.exists(phosphomimetic_e_folder_30):
                continue
            if not os.path.exists(phosphomimetic_d_folder_30):
                continue
            if not os.path.exists(wild_type_folder_50):
                continue
            if not os.path.exists(phosphomimetic_d_folder_50):
                continue
            if not os.path.exists(phosphomimetic_e_folder_50):
                continue

            for file in os.listdir(phosphomimetic_d_folder):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    phosphomimetic_d_pdb_full = os.path.join(phosphomimetic_d_folder, file)

            for file in os.listdir(wild_type_folder):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    wild_type_pdb_full = os.path.join(wild_type_folder, file)

            for file in os.listdir(phosphomimetic_e_folder):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    phosphomimetic_e_pdb_full = os.path.join(phosphomimetic_e_folder, file)

            for file in os.listdir(wild_type_folder_30):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    wild_type_pdb_60 = os.path.join(wild_type_folder_30, file)

            for file in os.listdir(phosphomimetic_d_folder_30):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    phosphomimetic_d_pdb_60 = os.path.join(phosphomimetic_d_folder_30, file)

            for file in os.listdir(phosphomimetic_e_folder_30):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    phosphomimetic_e_pdb_60 = os.path.join(phosphomimetic_e_folder_30, file)

            for file in os.listdir(wild_type_folder_50):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    wild_type_pdb_100 = os.path.join(wild_type_folder_50, file)

            for file in os.listdir(phosphomimetic_d_folder_50):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    phosphomimetic_d_pdb_100 = os.path.join(phosphomimetic_d_folder_50, file)

            for file in os.listdir(phosphomimetic_e_folder_50):
                if "_relaxed_rank_001_" in file and file.endswith(".pdb") and "sorted" not in file:
                    phosphomimetic_e_pdb_100 = os.path.join(phosphomimetic_e_folder_50, file)

            if os.path.isfile(wild_type_pdb_full) and os.path.isfile(phosphomimetic_d_pdb_full) and os.path.isfile(
                    phosphomimetic_e_pdb_full):
                sequence = extract_residues_from_PDB(wild_type_pdb_full)
                site = int(site)
                d = int(len(sequence))
                a = max(site - 5, 1)
                b = min(site + 5, d)
                e = max(site - 30, 1)
                f = min(site + 30, d)
                contigs_full = f"A1-{site - 1}/A{site + 1}-{d}"
                contigs_30 = f"A{a}-{site - 1}/A{site + 1}-{b}"
                contigs_50 = f"A{e}-{site - 1}/A{site + 1}-{f}"

                rmsd_d_full_contigs = CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_d_pdb_full)
                rmsd_e_full_contigs = CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_e_pdb_full)
                rmsd_d_30 = CompareTwoPDBs(contigs_30, wild_type_pdb_full, phosphomimetic_d_pdb_full)
                rmsd_e_30 = CompareTwoPDBs(contigs_30, wild_type_pdb_full, phosphomimetic_e_pdb_full)

                rmsd_d_full_contigs = round(rmsd_d_full_contigs, 3)
                rmsd_e_full_contigs = round(rmsd_e_full_contigs, 3)
                rmsd_d_30 = round(rmsd_d_30, 3)
                rmsd_e_30 = round(rmsd_e_30, 3)
                try:
                    rmsd_wt_full_vs_d_full = round(
                        CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_d_pdb_full), 3)
                except:
                    rmsd_wt_full_vs_d_full = 0

                try:
                    rmsd_wt_full_vs_e_full = round(
                        CompareTwoPDBs(contigs_full, wild_type_pdb_full, phosphomimetic_e_pdb_full), 3)
                except:
                    rmsd_wt_full_vs_e_full = 0

                try:
                    rmsd_wt_full_vs_wt_100 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_full, wild_type_pdb_100), 3)
                except:
                    rmsd_wt_full_vs_wt_100 = 0

                try:
                    rmsd_wt_full_vs_wt_60 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_full, wild_type_pdb_60), 3)
                except:
                    rmsd_wt_full_vs_wt_60 = 0

                try:
                    rmsd_d_full_vs_d_100 = round(
                        CompareTwoPDBs(contigs_30, phosphomimetic_d_pdb_full, phosphomimetic_d_pdb_100), 3)
                except:
                    rmsd_d_full_vs_d_100 = 0

                try:
                    rmsd_d_full_vs_d_60 = round(
                        CompareTwoPDBs(contigs_30, phosphomimetic_d_pdb_full, phosphomimetic_d_pdb_60), 3)
                except:
                    rmsd_d_full_vs_d_60 = 0

                try:
                    rmsd_e_full_vs_e_100 = round(
                        CompareTwoPDBs(contigs_30, phosphomimetic_e_pdb_full, phosphomimetic_e_pdb_100), 3)
                except:
                    rmsd_e_full_vs_e_100 = 0

                try:
                    rmsd_e_full_vs_e_60 = round(
                        CompareTwoPDBs(contigs_30, phosphomimetic_e_pdb_full, phosphomimetic_e_pdb_60), 3)
                except:
                    rmsd_e_full_vs_e_60 = 0

                try:
                    rmsd_wt_100_vs_d_100 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_100, phosphomimetic_d_pdb_100),
                                                 3)
                except:
                    rmsd_wt_100_vs_d_100 = 0

                try:
                    rmsd_wt_100_vs_e_100 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_100, phosphomimetic_e_pdb_100),
                                                 3)
                except:
                    rmsd_wt_100_vs_e_100 = 0

                try:
                    rmsd_wt_60_vs_d_60 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_60, phosphomimetic_d_pdb_60), 3)
                except:
                    rmsd_wt_60_vs_d_60 = 0

                try:
                    rmsd_wt_60_vs_e_60 = round(CompareTwoPDBs(contigs_30, wild_type_pdb_60, phosphomimetic_e_pdb_60), 3)
                except:
                    rmsd_wt_60_vs_e_60 = 0

                stability_wt = 0
                stability_d = 0
                stability_e = 0
                diff_score_D = 0
                diff_score_E = 0

                gene_data.extend([
                    rmsd_wt_full_vs_d_full, rmsd_wt_full_vs_e_full, rmsd_wt_full_vs_wt_100, rmsd_wt_full_vs_wt_60,
                    rmsd_d_full_vs_d_100, rmsd_d_full_vs_d_60, rmsd_e_full_vs_e_100, rmsd_e_full_vs_e_60,
                    rmsd_wt_100_vs_d_100, rmsd_wt_100_vs_e_100, rmsd_wt_60_vs_d_60, rmsd_wt_60_vs_e_60
                ])
                f_temp.write('\t'.join(map(str, gene_data)) + '\n')
                results.append(gene_data)
                print(f"{gene_name} done")
            else:
                print(f"Missing PDB files for gene: {gene_name}")

    os.rename(temp_file, output_file)
    # Sorting results according to RMSD D
    

    results_sorted_d = sorted(results, key=lambda x: sort_key_with_nan(x[-2]), reverse=True)
    write_sorted_results(output_file_sorted_d, header, results_sorted_d, use_pretty_table=True)

    # Sorting results according to RMSD E
    results_sorted_e = sorted(results, key=lambda x: sort_key_with_nan(x[-1]), reverse=True)
    write_sorted_results(output_file_sorted_e, header, results_sorted_e, use_pretty_table=True)

    results_sorted_e_diff = sorted(results, key=lambda x: sort_key_with_nan(x[-9]), reverse=True)
    write_sorted_results(output_file_sorted_e_diff, header, results_sorted_e_diff, use_pretty_table=True)

    results_sorted_d_diff = sorted(results, key=lambda x: sort_key_with_nan(x[-10]), reverse=True)
    write_sorted_results(output_file_sorted_d_diff, header, results_sorted_d_diff, use_pretty_table=True)

    # Writing results to text file using PrettyTable
    table = PrettyTable()
    table.field_names = header
    for result in results:
        table.add_row(result)

    # Remove borders and header repetition
    table.border = False
    table.header = True

    with open(output_txt, 'w') as f_txt:
        f_txt.write(table.get_string())

    # Writing results to CSV file
    with open(output_csv, 'w', newline='') as f_csv:
        writer = csv.writer(f_csv)
        writer.writerow(header)
        writer.writerows(results)

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
