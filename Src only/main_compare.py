from prettytable import PrettyTable
from confidence_score import find_pdb_file, extract_confidence_score
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from testbourré import calculate_stability_pdb
import csv
import os

# Define directories and input/output paths
input_file = '/mnt/c/Users/crisd/Desktop/Src/Src_in_vitro_tested.txt'
base_pdb_dir = '/mnt/c/Users/crisd/Desktop/Src/pdb'
output_file = '/mnt/c/Users/crisd/Desktop/Src/phosphorylation_analysis_results.txt'  # Path to save analysis results

# Définition des structures secondaires pour la comparaison
secondary_structures = [
    ('H', 6, 'alpha'), ('E', 6, 'beta'), ('B', 6, 'iso_b'),
    ('G', 6, 'alpha3'), ('I', 6, 'alphaI'), ('T', 6, 'hydturn')
]

def main_analysis(input_file, base_pdb_dir, output_file):
    # Charger les données de phosphorylation
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        phosphorylation_data = [row for row in reader]
    
    # Initialiser PrettyTable avec les nouvelles colonnes
    table = PrettyTable()
    table.field_names = [
        "Gene", "Site", "Conf_WT", "Conf_D", "Conf_E", "Diff_WT-D", "Diff_WT-E",
        "Diff_Score_WT-D", "Diff_Score_WT-E", "Diff_Alpha", "Diff_Beta", 
        "Diff_Iso_B", "Diff_Alpha3", "Diff_AlphaI", "Diff_Hydturn"
    ]
    
    # Stocker les résultats pour sauvegarde
    results = []
    
    for entry in phosphorylation_data:
        gene = entry['Gene']
        site = int(entry['Res'])
        
        # Rechercher et charger les fichiers PDB pour WT et phosphomimetic
        pdb_wt = find_pdb_file(gene, site, "", base_pdb_dir)
        pdb_D = find_pdb_file(gene, site, "D", base_pdb_dir)
        pdb_E = find_pdb_file(gene, site, "E", base_pdb_dir)
        
        # Calculer les scores de confiance
# Calculer les scores de confiance avec vérification de None avant d'arrondir
        conf_wt = round(extract_confidence_score(pdb_wt, site), 1) if pdb_wt and extract_confidence_score(pdb_wt, site) is not None else None
        conf_d = round(extract_confidence_score(pdb_D, site), 1) if pdb_D and extract_confidence_score(pdb_D, site) is not None else None
        conf_e = round(extract_confidence_score(pdb_E, site), 1) if pdb_E and extract_confidence_score(pdb_E, site) is not None else None


        # Différences de score de confiance
        diff_wt_d = round(conf_wt - conf_d, 1) if conf_wt and conf_d else None
        diff_wt_e = round(conf_wt - conf_e, 1) if conf_wt and conf_e else None

        if pdb_wt and pdb_D:
            try:
                wt_stability = calculate_stability_pdb(pdb_wt, [{"sub_site": site, "sub_letter": "Y"}], secondary_structures)
            except:
                wt_stability = None
            
            try:
                d_stability = calculate_stability_pdb(pdb_D, [{"sub_site": site, "sub_letter": "Y"}], secondary_structures)
            except:
                d_stability = None

            if wt_stability != None and d_stability != None:
                try:
                    stability_score_wt = wt_stability[0]['stability_score']
                    stability_score_d = d_stability[0]['stability_score']
                    diff_score_wt_d = stability_score_wt - stability_score_d
                except:
                    diff_score_wt_d = None
                
                try:
                    diff_alpha = wt_stability[0]['alpha'] - d_stability[0]['alpha']
                except:
                    diff_alpha = None
                
                try:
                    diff_beta = wt_stability[0]['beta'] - d_stability[0]['beta']
                except:
                    diff_beta = None
                
                try:
                    diff_iso_b = wt_stability[0]['iso_b'] - d_stability[0]['iso_b']
                except:
                    diff_iso_b = None
                
                try:
                    diff_alpha3 = wt_stability[0]['alpha3'] - d_stability[0]['alpha3']
                except:
                    diff_alpha3 = None
                
                try:
                    diff_alphaI = wt_stability[0]['alphaI'] - d_stability[0]['alphaI']
                except:
                    diff_alphaI = None
                
                try:
                    diff_hydturn = wt_stability[0]['hydturn'] - d_stability[0]['hydturn']
                except:
                    diff_hydturn = None
            else:
                diff_score_wt_d = diff_alpha = diff_beta = diff_iso_b = diff_alpha3 = diff_alphaI = diff_hydturn = None
        else:
            diff_score_wt_d = diff_alpha = diff_beta = diff_iso_b = diff_alpha3 = diff_alphaI = diff_hydturn = None

        # Stability scores and differences for each secondary structure for E
        if pdb_wt and pdb_E:
            try:
                e_stability = calculate_stability_pdb(pdb_E, [{"sub_site": site, "sub_letter": "Y"}], secondary_structures)
            except:
                e_stability = None

            if wt_stability != None and e_stability != None:
                try:
                    stability_score_e = e_stability[0]['stability_score']
                    diff_score_wt_e = stability_score_wt - stability_score_e
                except:
                    diff_score_wt_e = None

                try:
                    diff_alpha_e = wt_stability[0]['alpha'] - e_stability[0]['alpha']
                except:
                    diff_alpha_e = None

                try:
                    diff_beta_e = wt_stability[0]['beta'] - e_stability[0]['beta']
                except:
                    diff_beta_e = None

                try:
                    diff_iso_b_e = wt_stability[0]['iso_b'] - e_stability[0]['iso_b']
                except:
                    diff_iso_b_e = None

                try:
                    diff_alpha3_e = wt_stability[0]['alpha3'] - e_stability[0]['alpha3']
                except:
                    diff_alpha3_e = None

                try:
                    diff_alphaI_e = wt_stability[0]['alphaI'] - e_stability[0]['alphaI']
                except:
                    diff_alphaI_e = None

                try:
                    diff_hydturn_e = wt_stability[0]['hydturn'] - e_stability[0]['hydturn']
                except:
                    diff_hydturn_e = None
            else:
                diff_score_wt_e = diff_alpha_e = diff_beta_e = diff_iso_b_e = diff_alpha3_e = diff_alphaI_e = diff_hydturn_e = None
        else:
            diff_score_wt_e = diff_alpha_e = diff_beta_e = diff_iso_b_e = diff_alpha3_e = diff_alphaI_e = diff_hydturn_e = None

        # Ajouter les résultats dans le tableau PrettyTable
        table.add_row([
            gene, site, conf_wt, conf_d, conf_e, diff_wt_d, diff_wt_e,
            diff_score_wt_d, diff_score_wt_e, diff_alpha, diff_beta, 
            diff_iso_b, diff_alpha3, diff_alphaI, diff_hydturn
        ])
        
        # Ajouter aux résultats pour sauvegarde
        results.append([
            gene, site, conf_wt, conf_d, conf_e, diff_wt_d, diff_wt_e,
            diff_score_wt_d, diff_score_wt_e, diff_alpha, diff_beta, 
            diff_iso_b, diff_alpha3, diff_alphaI, diff_hydturn
        ])
    
    with open(output_file, 'w') as f_out:

        f_out.write(table.get_string())  # Écrire PrettyTable sous forme de texte

        # Écrire les résultats dans un format tabulaire
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow([
            "Gene", "Site", "Conf_WT", "Conf_D", "Conf_E", "Diff_WT-D", "Diff_WT-E",
            "Diff_Score_WT-D", "Diff_Score_WT-E", "Diff_Alpha", "Diff_Beta", 
            "Diff_Iso_B", "Diff_Alpha3", "Diff_AlphaI", "Diff_Hydturn"
        ])
        writer.writerows(results)

main_analysis(input_file, base_pdb_dir, output_file)