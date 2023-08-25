from  import process_gene_site
from no_id_given.graph_no_id_given import print_outside_file

# This branch allows the analyse of the phosphorylation sites given one site_list where the id of the protein is not given.
# Before any use, you should adapt the function create_gene_dict in order to process the information on your site_list correctly, and maybe the check_phosphorylation depending on what information do you have on your site_list
# Replace these paths with the actual paths for your system

site_list_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/yeast/PKA/site_list.txt"  # Replace with the actual path to site_list.txt
cif_directory = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/yeast_database/UP000002311_559292_YEAST_v4"  # Replace with the actual path to the directory containing CIF files
output_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/yeast/PKA/results.txt"  # Replace with the desired output file path for all datas
output_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/yeast/PKA/results"  # Replace with the desired ouptut folder for the stats
database_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/yeast/Cdk1/DatabaseToCompare"  # Replace with the folder containing the databases with which you want to compare. The function search_phosphorylation_info may need to be adapt
output_file_path = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/yeast/PKA/stats_output.txt"  # Replace with the desired output file path

secondary_structure_distance=8
count=0
secondary_structures=[('H',secondary_structure_distance,'alpha',count),('E',secondary_structure_distance,'beta',count),('B',secondary_structure_distance,'iso_B',count),('G',secondary_structure_distance,'alpha3',count),('I',secondary_structure_distance,'alphaI',count),('T',secondary_structure_distance,'hydturn',count)]
to_rank=[('Score',15,False),('Confidence',15,False),('NextAA',15,True),('alpha',15,False),('beta',15,False),('iso_B',15,False),('alpha3',15,False),('alphaI',15,False),('hydturn',15,False)]

process_gene_site(site_list_file,cif_directory,output_file,database_folder,secondary_structures)
print_outside_file(output_file_path,output_file,to_rank,secondary_structures,output_folder)
