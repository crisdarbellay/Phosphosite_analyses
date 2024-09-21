import os
import csv

def check_missing_predictions(selection_file, fastas_folder, missing_predictions_file):
    seen_entries = set()
    missing_entries = []

    # Open the selection file and process it
    with open(selection_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        header = next(reader)  # Read and skip the header

        for row in reader:
            gene = row[1]
            site = row[2]
            entry = (gene, site)
            fasta_filename = f"{gene}-D{site}.fasta"
            fasta_file_path = os.path.join(fastas_folder, fasta_filename)

            # Check if we have seen this entry before
            if entry not in seen_entries:
                seen_entries.add(entry)
                # Check if the file exists
                if not os.path.isfile(fasta_file_path):
                    missing_entries.append(row)

    # Write the missing entries to the output file
    with open(missing_predictions_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['gene', 'site'] + header[2:])  # Write header with 'gene', 'site' and other columns
        writer.writerows(missing_entries)

# Example usage
selection_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/2_algorithm/selection.txt"
fastas_folder = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/2_algorithm/phosphosites/fastas"
missing_predictions_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/2_algorithm/missing_predictions.txt"
check_missing_predictions(selection_file, fastas_folder, missing_predictions_file)

print("Missing predictions check completed.")
