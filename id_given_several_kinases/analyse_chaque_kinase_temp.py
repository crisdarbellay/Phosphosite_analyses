import pandas as pd
from collections import defaultdict

def load_kinase_data(kinase_file):
    kinase_data = pd.read_csv(kinase_file, delimiter='\t')
    kinase_data = kinase_data[(kinase_data['KIN_ORGANISM'] == 'human') & (kinase_data['SUB_ORGANISM'] == 'human')]
    kinase_dict = defaultdict(list)
    for _, row in kinase_data.iterrows():
        kinase_dict[(row['SUB_GENE'], str(row['SUB_MOD_RSD'][1:]))].append(row['KINASE'])
    return kinase_dict

def process_gene_site_data(gene_site_file, kinase_data, output_file):
    gene_site_data = pd.read_csv(gene_site_file, delim_whitespace=True)
    gene_site_data['Kinases'] = gene_site_data.apply(lambda row: ", ".join(kinase_data.get((row['Gene'], str(row['site'])), [])), axis=1)
    
    kinase_counts = defaultdict(lambda: defaultdict(int))
    kinase_pairs = defaultdict(list)
    
    for _, row in gene_site_data.iterrows():
        gene = row['Gene']
        site = row['site']
        kinases = kinase_data.get((gene, str(site)), [])
        for kinase in kinases:
            kinase_counts[kinase][(gene, site)] += 1
            kinase_pairs[kinase].append((gene, site))
    
    with open(output_file, 'w', newline='') as kinfile:
        kinfile.write("Kinase\tgenes_tot\tGene\tsite\n")
        for kinase, pairs in kinase_pairs.items():
            genes_tot = len(pairs)
            for gene, site in pairs:
                kinfile.write(f"{kinase}\t{genes_tot}\t{gene}\t{site}\n")

# Exemple d'utilisation
kinase_file = r"/mnt/c/Users/crisd/Downloads/Kinase_Substrate_Dataset/Kinase_Substrate_Dataset"
gene_site_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/to_test/selection_full.txt"
output_file = r"/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/to_test/kinases_paires.txt"

kinase_data = load_kinase_data(kinase_file)
process_gene_site_data(gene_site_file, kinase_data, output_file)

print("Processing completed.")
