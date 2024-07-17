import os

def extract_info_from_file(file_path):
    phosphosites = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        ranking_score_started = False
        ranking_alpha_started = False
        for line in lines:
            if line.startswith('Ranking Score:'):
                ranking_score_started = True
                ranking_alpha_started = False
            elif line.startswith('Ranking Alpha:'):
                ranking_alpha_started = True
                ranking_score_started = False
            elif ranking_score_started and not ranking_alpha_started and line.strip():
                try:
                    data = line.strip().split('\t')
                    gene = data[1]
                    site = data[2]
                except:
                    print(line)
                # Vérifier si la paire gene,site est déjà dans le dictionnaire
                if (gene, site) not in phosphosites and site.isdigit():
                    phosphosites[(gene, site)] = data[3:10]
    return phosphosites

def get_top_phosphosites(folder_path):
    top_phosphosites = {}
    for subdir, _, files in os.walk(folder_path):
        kinase_name = os.path.basename(subdir).replace('_human_data', '')
        if kinase_name:
            top_phosphosites[kinase_name] = []
            for file in files:
                if file.endswith('.txt'):
                    file_path = os.path.join(subdir, file)
                    phosphosites = extract_info_from_file(file_path)
                    for pair, info in phosphosites.items():
                        top_phosphosites[kinase_name].append([pair[0], pair[1]] + info)
    return top_phosphosites

def filter_phosphosites(top_phosphosites):
    filtered_phosphosites = {}
    for kinase, phosphosites in top_phosphosites.items():
        unique_phosphosites = set()
        filtered_phosphosites[kinase] = []
        for phosphosite in phosphosites:
            gene = phosphosite[0]
            site = phosphosite[1]
            if (gene, site) not in unique_phosphosites:
                unique_phosphosites.add((gene, site))
                filtered_phosphosites[kinase].append(phosphosite)
    return filtered_phosphosites

def write_to_file(filtered_phosphosites):
    with open('/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results.txt', 'w') as file:
        file.write("Gene\tsite\tScore\tAFConf\tNextAA\tin_vivo\tin_vitro\tid_alphafold\tsite_aa\n")
        for kinase, phosphosites in filtered_phosphosites.items():
            file.write(f"Kinase {kinase} :\n")
            for phosphosite in phosphosites:
                file.write('\t'.join(phosphosite) + '\n')

# Dossier contenant les données des kinases
folder_path = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/paper/kinase_2_substrates/human/results"

# Obtention des meilleurs phosphosites pour chaque kinase
top_phosphosites = get_top_phosphosites(folder_path)


# Écrire les informations filtrées dans un fichier
write_to_file(top_phosphosites)
