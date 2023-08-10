from utility import process_gene_site
from graph import generate_3d_surface_graphs

                
# Replace these paths with the actual paths for your system
site_list_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/Cdk1/site_list.txt"  # Replace with the actual path to site_list.txt
cif_directory = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/yeast_database/UP000002311_559292_YEAST_v4"  # Replace with the actual path to the directory containing PDB and CIF files
output_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/results.txt"  # Replace with the desired output file path
database_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/Cdk1/DatabaseToCompare"
output_folder = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/Cdk1/Graph"
process_gene_site(site_list_file,cif_directory,output_file,database_folder)
generate_3d_surface_graphs(output_file,output_folder)