from utility import process_gene_site
from graph import print_outside_file

                
# Replace these paths with the actual paths for your system
site_list_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/PKA/site_list.txt"  # Replace with the actual path to site_list.txt
cif_directory = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/yeast_database/UP000002311_559292_YEAST_v4"  # Replace with the actual path to the directory containing PDB and CIF files
output_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/PKA/results.txt"  # Replace with the desired output file path
database_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/Cdk1/DatabaseToCompare"
output_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/PKA/results"
output_file_path = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/PKA/results/stats_output.txt"  # Replace with the desired output file path

alpha_H_distance=10
beta_E_distance=10
count=0
secondary_structures=[('H',alpha_H_distance,'alpha',count),('E',beta_E_distance,'beta',count),('B',beta_E_distance,'iso_B',count),('G',alpha_H_distance,'alpha3',count),('I',alpha_H_distance,'alphaI',count),('T',alpha_H_distance,'hydrogene_turn',count)]
to_rank=[('Score',10,False),('Confidence',10,False),('NextAA',10,True),('alpha',10,False),('beta',10,False),('iso_B',10,False),('alpha3',10,False),('alphaI',10,False),('hydrogene_turn',10,False)]

process_gene_site(site_list_file,cif_directory,output_file,database_folder,secondary_structures)
print_outside_file(output_file_path,output_file,to_rank,secondary_structures)
