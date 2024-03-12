import os
from graph_several_kinases import print_outside_file
from id_given_several_kinases_analysis import create_kinases_datas

# This branch allows the analyse of the phosphorylation sites given one site_list, containing several kinases and several substrates, and where the id of the protein is given.
# Before any use, you should adapt the function create_gene_dict in order to process the information on your site_list correctly, and maybe the check_phosphorylation depending on what information do you have on your site_list
# Replace these paths with the actual paths for your system


human_cif_directory = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/human_database"                #path to your human cif database folder (download on AlphaFold website)
mouse_cif_directory = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/mouse_database"                #path to your mouse cif database folder         ""
rat_cif_directory = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/rat_database"                    #path to your rat cif database folder           ""
output_directory = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/kinase_2_substrates"                          #path to your results folder
file_path = r"/mnt/c/Users/crisd/Downloads/Kinase_Substrate_Dataset/Kinase_Substrate_Dataset"    #path to your site_list
input_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/kinase_2_substrates/human"                        #path to the input folder for which you want to generate datas (human,mouse,rat)
output_root_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results"          #path to the output folder for your stats

sec_struct_distance = 8     #distance of detection for the secondary structures around your sites

count = 0
secondary_structures = [('H', sec_struct_distance, 'alpha', count), ('E', sec_struct_distance, 'beta', count), ('B', sec_struct_distance, 'iso_B', count), ('G', sec_struct_distance, 'alpha3', count), ('I', sec_struct_distance, 'alphaI', count), ('T', sec_struct_distance, 'hydturn', count)]
to_rank = [('Score', 15, False), ('AFConf', 15, False), ('NextAA', 15, True), ('alpha', 15, False), ('beta', 15, False),('hydturn',15,False)]

create_kinases_datas(file_path,secondary_structures,human_cif_directory,mouse_cif_directory,rat_cif_directory,output_directory)

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