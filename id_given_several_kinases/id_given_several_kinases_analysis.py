import pandas as pd
import os
import os
from utility_id_given import calculate_stability_cif


def convert_code_to_filename(code):
    return f"AF-{code}-F1-model_v4.cif"

def create_kinase_dictionary(file_path):
    kinase_data = {}

    # Lire le fichier CSV en utilisant pandas
    df = pd.read_csv(file_path, sep='\t', encoding='iso-8859-1')

    for _, row in df.iterrows():
        kinase = row['KINASE']+'_'+row['KIN_ORGANISM']
        kinase_animal=row['KIN_ORGANISM']
        substrate_gene = row['SUB_GENE']
        acc_id = row['KIN_ACC_ID']
        substrate = row['SUBSTRATE']
        substrate_gene_id = row['SUB_GENE_ID']
        substrate_mod_rsd = row['SUB_MOD_RSD']
        site_group_id = row['SITE_GRP_ID']
        site_aa = row['SITE_+/-7_AA']
        domain = row['DOMAIN']
        in_vivo = 'X' in row['IN_VIVO_RXN']
        in_vitro = 'X' in row['IN_VITRO_RXN']
        sub_acc_id=row['SUB_ACC_ID']
        sub_organisme=row['SUB_ORGANISM']
        

        if kinase not in kinase_data:
            kinase_data[kinase] = {'kinase_animal':kinase_animal}

        if substrate_gene not in kinase_data[kinase]:
            kinase_data[kinase][substrate_gene] = []

        kinase_data[kinase][substrate_gene].append({
            'acc_id': acc_id,
            'substrate': substrate,
            'substrate_gene':substrate_gene,
            'substrate_gene_id': substrate_gene_id,
            'sub_acc_id':sub_acc_id,
            'sub_site': substrate_mod_rsd[1:],
            'sub_letter':substrate_mod_rsd[0],
            'site_group_id': site_group_id,
            'kinase_animal':kinase_animal,
            'site_aa': site_aa,
            'domain': domain,
            'in_vivo': in_vivo,
            'in_vitro': in_vitro,
            'filename':convert_code_to_filename(sub_acc_id),
            'sub_organism':sub_organisme
            })

    return kinase_data

def create_kinases_datas(file_path,secondary_structures,human_cif_directory,mouse_cif_directory,rat_cif_directory,output_directory):

    kinase_dictionary = create_kinase_dictionary(file_path)

    i=0

    for kinase, gene_data in kinase_dictionary.items():
        results = []
        i=i+1
        print(f"{kinase}   {i}")
    # Create a folder for each kinase
    
        output_directory_=output_directory+f"/{gene_data['kinase_animal']}"
        output_file = os.path.join(output_directory_, f"{kinase}_results.txt")


            # Define name of files 
        if not os.path.exists(output_directory_):
            os.makedirs(output_directory+f"/{gene_data['kinase_animal']}")   
            output_file = os.path.join(output_directory_, f"{kinase}_results.txt")
        
    # Separate group in_vivo and in_vitro
        for index, (gene, sites_data) in enumerate(gene_data.items()):
        
            if index ==0:    
                continue  
            organism=sites_data[0]['sub_organism']
            if organism=='human':
                cif_file = human_cif_directory + '/' + sites_data[0]['filename'] + '.gz'
            elif organism=='mouse':
                cif_file = mouse_cif_directory + '/' + sites_data[0]['filename'] + '.gz'
            elif organism=='rat':
                cif_file = rat_cif_directory + '/' + sites_data[0]['filename'] + '.gz'   
            else:
                continue 
            try:
                stability_scores = calculate_stability_cif(cif_file, sites_data, secondary_structures)
                results.append({'gene_name': gene, 'stability_scores': stability_scores, 'sub_acc_id':sites_data[0]['sub_acc_id']})
            except:
                continue

        sorted_results = sorted(results, key=lambda x: x['gene_name'])


        with open(output_file, 'w') as out_file:
            header = "Gene\tRes\tScore\talpha\tbeta\tiso_B\talpha3\talphaI\thydturn\tAFConf\tNextAA\torga\tsite_aa\tsub_id\tin_vivo\tin_vitro\n"
            out_file.write(header)

            for result in sorted_results:
                gene_name = result['gene_name']
                stability_scores = result['stability_scores']
                sub_acc_id=result['sub_acc_id']

                for stats in stability_scores:
                    formatted_values = [
                        gene_name,
                        stats['residue_number'],
                        stats['stability_score'],
                        stats['alpha'],
                        stats['beta'],
                        stats['iso_B'],
                        stats['alpha3'],
                        stats['alphaI'],
                        stats['hydturn'],
                        stats['average_confidence_score'],
                        stats['chain_dist_secondary_struct'],
                        stats['sub_organism'],
                        stats['site_aa'],
                        sub_acc_id,
                        "X" if stats['in_vivo'] else "-",
                        "X" if stats['in_vitro'] else "-"
                ]

            # Write values in output file
                    formatted_line = "\t".join(map(str, formatted_values)) + "\n"
                    out_file.write(formatted_line)


  