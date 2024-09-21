def filter_txt(input_file, output_file):
    """
    This function filters the rows of a TXT file based on specific criteria and writes the filtered rows to a new file.
    
    Criteria:
    - Keep the row if at least one of the values in the columns `wt30-D30_5` or `wt30-E30_5` (line[-4], line[-3]) is greater than or equal to 2.
    - OR, keep the row if at least one of the values in the columns `wt30-E30_5` or `wt30-D30_5` (line[-1], line[-2]) is greater than or equal to 1.
    
    :param input_file: Path to the input TXT file.
    :param output_file: Path to the output filtered TXT file.
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = infile.readlines()

        header = lines[0]  # Read the header row
        outfile.write(header)  # Write the header to the output file

        for line in lines[1:]:
            # Split the line by tabs
            columns = line.strip().split('\t')
            
            # Convert the relevant columns to floats and check the criteria
            if (float(columns[-4]) >= 2 or float(columns[-3]) >= 2) or (float(columns[-1]) >= 1 or float(columns[-2]) >= 1):
                outfile.write(line)

# Exemple d'utilisation
input_file = '/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_rmsd.txt'  # Remplacez par le chemin vers votre fichier d'entrée
output_file = '/mnt/c/Users/crisd/Desktop/Phosphoswitch/datas/6_angstrom/4_algorithm/results/results_rmsd_filtered.txt'  # Remplacez par le chemin où vous souhaitez enregistrer le fichier filtré
filter_txt(input_file, output_file)
