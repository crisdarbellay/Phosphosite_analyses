import subprocess
import os
import tempfile

def count_alpha_helices(pdb_path):
    temp_dssp_file = tempfile.NamedTemporaryFile(delete=False)
    temp_dssp_file_path = temp_dssp_file.name
    temp_dssp_file.close()

    wsl_dssp_command = f"wsl mkdssp {pdb_path} {temp_dssp_file_path}"
    subprocess.call(wsl_dssp_command, shell=True)

    alpha_helix_count = 0

    with open(temp_dssp_file_path, 'r') as dssp_file:
        parsing_helix = False
        for line in dssp_file:
            if line.startswith("  #  RESIDUE AA STRUCTURE"):
                parsing_helix = True
            elif parsing_helix and line.startswith("  #"):
                parsing_helix = False
            elif parsing_helix:
                helix_info = line.split()
                helix_structure = helix_info[2]
                if helix_structure == "H":
                    alpha_helix_count += 1

    os.remove(temp_dssp_file_path)

    return alpha_helix_count

pdb_path = r"C:\Users\crisd\Desktop\Projet\NTH3.pdb"
num_alpha_helices = count_alpha_helices(pdb_path)
print(f"Le nombre d'hélices alpha dans NTH3 est : {num_alpha_helices}")
