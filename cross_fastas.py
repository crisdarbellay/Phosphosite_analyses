import os
import shutil

def copy_missing_fastas(human_folder, predicted_folder, to_predict_folder,too_long_file):
    if not os.path.exists(to_predict_folder):
        os.makedirs(to_predict_folder)

    human_fastas = set(os.listdir(human_folder))
    predicted_fastas = set(os.listdir(predicted_folder))

    for fasta in human_fastas:
        if fasta not in predicted_fastas:
            src_path = os.path.join(human_folder, fasta)
            dst_path = os.path.join(to_predict_folder, fasta)
            shutil.copy(src_path, dst_path)
            print(f"Copied {fasta} to {to_predict_folder}")
    with open(too_long_file, 'r') as file:
        lines = file.readlines()
        header = lines[0]
        lines=lines[1:]
        for line in lines:
            gene = line.split()[1]
            site = line.split()[2]
            path1 = f"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/all_fastas_to_test/{gene}-{site}.fasta"
            path2 =f"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/all_fastas_to_test/{gene}-D{site}.fasta"
            path3 = f"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/all_fastas_to_test/{gene}-E{site}.fasta"
            shutil.copy(path1, dst_path)
            shutil.copy(path2, dst_path)
            shutil.copy(path3, dst_path)


# Example usage
human_folder = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/all_fastas_to_test"
predicted_folder = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/8_angstrom/2_algorithm/phosphosites/fastas-predicted"
to_predict_folder = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/to_predict"
too_long_file = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/2_algorithm/phosphosites/too_long.txt"

copy_missing_fastas(human_folder, predicted_folder, to_predict_folder,too_long_file)

print("Copying completed.")
