import os
import re

def process_pdb_file(pdb_file, site, output_folder):
    output_path = os.path.join(output_folder, f'{os.path.basename(pdb_file).split(".")[0]}.pdb')
    decalage = max(0, site - 50)
    if decalage > 0:
        with open(pdb_file, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    parts = line.split()
                    res_num = int(parts[5])
                    new_res_num = res_num + decalage
                    parts[5] = str(new_res_num)
                    new_line = '\t'.join(parts)
                    outfile.write(new_line + '\n')
    
    print(f"Processed {pdb_file}, output saved to {output_path}")

def split_string(s):
    # Définir le modèle regex pour correspondre aux trois formes différentes
    pattern = re.compile(r'^(?P<protein>.+?)-(?:D|E)?(?P<site>.+)$')
    match = pattern.match(s)
    if match:
        protein = match.group('protein')
        try:
            site = int(match.group('site'))
        except:
            return 'H1-2',55
        return protein, site


output_folder = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-50-contigs-adapted/"
base_folder = f"/mnt/c/Users/crisd/Desktop/ProteinDesign/mass_test_1/pdb-ball-50"

for root, dirs, files in os.walk(base_folder):
    for file in files:
        if file.endswith(".pdb") and "_relaxed_rank_001_alphafold2_" in file and "sorted" not in file:
            protein_site = root.split("/")[-1]
            output_folder_temp = output_folder + protein_site
            protein,site =split_string(protein_site)
            file = root + '/' + file 
            process_pdb_file(file, site,output_folder_temp)
