import os
from graph import print_outside_file

# Define the input and output folders
input_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/results/human"
output_root_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/results/human/results"

alpha_H_distance = 10
beta_E_distance = 10
count = 0
secondary_structures = [('H', alpha_H_distance, 'alpha', count), ('E', beta_E_distance, 'beta', count), ('B', beta_E_distance, 'iso_B', count), ('G', alpha_H_distance, 'alpha3', count), ('I', alpha_H_distance, 'alphaI', count), ('T', alpha_H_distance, 'hydrogene_turn', count)]
to_rank = [('Score', 15, False), ('AFConf', 15, False), ('NextAA', 15, True), ('alpha', 15, False), ('beta', 15, False),('hydturn',15,False)]

def process_and_save(file_path, output_folder):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    output_file_path = os.path.join(output_folder, "output_file.txt")
    print_outside_file(output_file_path, file_path, to_rank, secondary_structures,output_folder,kinase="")
    
# Iterate through files in the input folder
for file_name in os.listdir(input_folder):
    if file_name.endswith(".txt"):
        input_file_path = os.path.join(input_folder, file_name)
        output_subfolder_name = file_name.replace("_results.txt", "_datas")
        output_subfolder_path = os.path.join(output_root_folder, output_subfolder_name)
        
        process_and_save(input_file_path, output_subfolder_path)
