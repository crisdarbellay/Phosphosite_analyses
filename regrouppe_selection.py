import os
import csv
from collections import defaultdict

def identify_and_remove_duplicates(input_folder, output_file, output_no_duplicates_file):
    entry_count = defaultdict(int)
    all_lines = []

    # Step 1: Identify and count duplicates
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".txt"):
            kinase_name = file_name.replace('_human_results.txt', '')
            input_file_path = os.path.join(input_folder, file_name)
            
            with open(input_file_path, 'r') as infile:
                reader = csv.reader(infile, delimiter='\t')
                header = next(reader)  # Skip header

                for row in reader:
                    gene = row[0]
                    site = row[1]
                    entry = (gene, site)
                    entry_count[entry] += 1
                    all_lines.append((kinase_name, row))

    # Step 2: Write the output file including duplicates
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['kinase'] + header)  # Write header

        for kinase_name, row in all_lines:
            gene = row[0]
            site = row[1]
            entry = (gene, site)
            writer.writerow([kinase_name] + row)

    # Step 3: Write the output file excluding duplicates
    with open(output_no_duplicates_file, 'w', newline='') as outfile_no_dup:
        writer = csv.writer(outfile_no_dup, delimiter='\t')
        writer.writerow(['kinase'] + header)  # Write header

        for kinase_name, row in all_lines:
            gene = row[0]
            site = row[1]
            entry = (gene, site)
            if entry_count[entry] == 1:  # Only write non-duplicate entries
                writer.writerow([kinase_name] + row)

# Example usage
input_folder = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/control/1_algorithm/all_kinases_results_filtered/human"
output_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/control/2_algorithm/selection.txt"
output_no_duplicates_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/control/2_algorithm/selection_no_doublon.txt"
identify_and_remove_duplicates(input_folder, output_file, output_no_duplicates_file)

print("Data processing completed.")
