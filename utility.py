import tempfile
import math
import subprocess
import os
import zipfile
import gzip
from Bio.PDB import PDBParser
import shutil
from Bio.PDB import MMCIFParser
import io
import matplotlib.pyplot as plt


def begin_by_ATOM(pdb_path):
    """
    Extract the lines starting with 'ATOM' from a PDB file.
    """
    atom_lines = []
    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith('ATOM'):
                atom_lines.append(line)
    return atom_lines

def calculate_distance(coord1, coord2):
    """
    Calculate the Euclidean distance between two 3D coordinates.
    """
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    distance = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    return distance

def extract_confidence_score(line):
    """
    Extract the confidence score from a line in the PDB file.
    """
    confidence_score = float(line.split()[14])
    return confidence_score



def search_phosphorylation_info(gene_dict, database_folder):
    base_databases = os.listdir(database_folder)
    for database_file in base_databases:
        database_path = os.path.join(database_folder, database_file)
        with open(database_path, 'r') as db_file:
            for line in db_file:
                for gene_name, sites in gene_dict.items():
                    if gene_name in line and 'phosphorylated' in line:
                        for site_number, site_info in enumerate(sites):
                            if str(site_info[0]) in line:
                                # Add a new tuple with the database name and True
                                gene_dict[gene_name][site_number] += ((database_file, True),)
                                break
    return gene_dict





def calculate_confidence_scores(atom_line, residue_key):
    """
    Calculate the confidence scores for a specific residue.
    """
    confidence_scores = {}
    for line in atom_line:           
        line_residue_number = int(line.split()[-2])            
        if line_residue_number == residue_key[1]:
            confidence_score = extract_confidence_score(line)
            confidence_scores[residue_key] = confidence_score
            break
    return confidence_scores


def calculate_stability_cif(cif_file_path, site):
    """
    Calculate the stability scores for residues at specific locations in a CIF file.
    """
    is_compressed = cif_file_path.endswith(".gz")

    # If the file is compressed, create a temporary decompressed file
    if is_compressed:
        temp_cif_file_path = tempfile.mkstemp(suffix=".cif")[1]
        with gzip.open(cif_file_path, 'rt') as gzipped_cif, open(temp_cif_file_path, 'w') as temp_cif:
            shutil.copyfileobj(gzipped_cif, temp_cif)
        cif_file_path = temp_cif_file_path

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", cif_file_path)
    model = structure[0]  # Assuming there is only one model in the structure

    stability_scores = []

    # Generate DSSP file
    temp_dssp_file_path = tempfile.mkstemp()[1]

    if os.name == 'nt':
        # Windows environment
        wsl_dssp_command = f"wsl mkdssp {cif_file_path} {temp_dssp_file_path}"
        subprocess.call(wsl_dssp_command, shell=True)
    else:
        # Non-Windows environment
        dssp_executable = "mkdssp"
        subprocess.call([dssp_executable, cif_file_path, temp_dssp_file_path,"--output-format=dssp"]) 

    # Store DSSP information in a dictionary
    dssp_info = {}
    with open(temp_dssp_file_path, 'r') as dssp_file:
        dssp_lines = dssp_file.readlines()
        for dssp_line in dssp_lines:
            line=dssp_line.split()
            if line[1]==line[0] and line[2]=='A':
                residue_number = int(line[1])
                res_dssp = line[4]
                dssp_info[residue_number] = res_dssp              

    helix_count = 0
    beta_count = 0
    stability_score = 0
    chain_distance = 100000
    chain_distance_temp = 10000
    confidence_scores_sum = 0
    confidence_scores_count = 0
    alpha_helix_distance = 10.0
    beta_fold_distance = 10.0

    residues = model.get_residues()
    for residu1 in residues:
        if residu1.get_id()[1] == site:
            for residu2 in model.get_residues():
                distance = calculate_distance(residu1['CA'].get_coord(), residu2['CA'].get_coord())
                res_dssp = dssp_info.get(residu2.get_id()[1], '')

                if distance <= alpha_helix_distance:                        
                    if res_dssp == 'H':
                        helix_count += 1
                if distance <= beta_fold_distance:
                    if res_dssp == 'E':
                        beta_count += 1
                if distance <= 10.0:
                    try:
                        confidence_score = calculate_confidence_scores(begin_by_ATOM(temp_cif_file_path), residu2.get_id())[residu2.get_id()]
                    except: 
                        confidence_score = None
                    if confidence_score is not None:
                        confidence_scores_sum += confidence_score
                        confidence_scores_count += 1                           
                if res_dssp == 'E' or res_dssp == 'H':
                    chain_distance_temp = abs(residu1.get_id()[1] - residu2.get_id()[1])
                    if chain_distance_temp < chain_distance:
                        chain_distance = chain_distance_temp

    stability_score = helix_count + beta_count
    average_confidence_score = round(confidence_scores_sum / confidence_scores_count, 1) if confidence_scores_count > 0 else None

    stability_scores.append({
            'residue_number': site,
            'stability_score': stability_score,
            'helix_count': helix_count,
            'beta_count': beta_count,
            'average_confidence_score': average_confidence_score,
            'chain_dist_secondary_struct': chain_distance
        })

    # Remove temporary DSSP file
    os.remove(temp_dssp_file_path)

    return stability_scores



def extract_protein_sequence_from_cif(cif_file):
    """
    Extract the protein sequence from a CIF file.
    
    """
    protein_sequence = ""
    amino_acid_code = None

    for line in cif_file:
        if line.startswith("ATOM"):
            parts = line.split()
            if len(parts) >= 6 and parts[3] == 'CA':
                amino_acid_code = parts[-1]
                protein_sequence += amino_acid_code

    return protein_sequence


def check_phosphorylation_context(sequence, target_position, context):
    """
    Check if the phosphorylation context around the target position in the sequence matches the given context.
    
    """
    
    start_position = target_position - 3
    end_position = target_position + 2
    
    target_sequence = sequence[start_position:end_position]
    
    return target_sequence == context


def create_gene_dict(site_list_file):
    gene_dict = {}
    with open(site_list_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            gene_name = parts[1]
            site_position = int(parts[2])
            phosphorylation_context_gross = parts[3]
            position_hashtag=phosphorylation_context_gross.index('#')
            phosphorylation_context_pre=phosphorylation_context_gross[position_hashtag-3:position_hashtag+3]
            phosphorylation_context=phosphorylation_context_pre.replace('#','')
            
            if gene_name not in gene_dict:
                gene_dict[gene_name] = []
            gene_dict[gene_name].append((site_position, phosphorylation_context))
    return gene_dict

def process_gene_site(site_list_file, cif_directory, output_file,database_folder):
    """
    Process the gene sites from the site list file, find corresponding CIF and PDB files,
    and calculate stability for each site. Write the results to the output file in a tab-separated format.
    """
    results = []
    base_databases = [file for file in os.listdir(database_folder) if file.endswith('.txt')]   
    # First, process CIF files and extract gene names
    gene_dict = create_gene_dict(site_list_file)  
    phosphorylation_info = search_phosphorylation_info(gene_dict, database_folder)
    
    for file_name in os.listdir(cif_directory):
        print(file_name)
        if file_name.endswith(".cif.gz"):
            cif_file_path = os.path.join(cif_directory, file_name)
            with gzip.open(cif_file_path, 'rt') as cif_file:
                key_residue=None
                for line in cif_file:
                    if "_ma_target_ref_db_details.gene_name" in line:
                        for key in gene_dict.keys():
                            if key in line:
                                key_residue=key
                if not key_residue==None:
                    cif_file.seek(0)

                    protein_sequence=extract_protein_sequence_from_cif(cif_file)
                    for site in gene_dict[key_residue]:
                        if check_phosphorylation_context(protein_sequence, site[0], site[1]):
                            stability_scores = calculate_stability_cif(cif_file_path, site[0])
                            results.append({
                                'gene_name': key_residue,
                                'stability_scores': stability_scores,
                                'phosphorylation_info': phosphorylation_info
                            })
                         
                            
    sorted_results = sorted(results, key=lambda x: x['gene_name'])

# Write sorted results to the output file
    with open(output_file, 'w') as out_file:    
        out_file.write("Gene\tRes\tScore\tHelix\tBeta\tConfidence\tNextAA")
        for database in base_databases:
         out_file.write(f"\t{database}")
        out_file.write("\n")
    
        for result in sorted_results:
            gene_name = result['gene_name']
            stability_scores = result['stability_scores']
            phosphorylation_info = result['phosphorylation_info']
        
            for stats in stability_scores:
                residue_number = stats['residue_number']
                stability_score = stats['stability_score']
                helix_count = stats['helix_count']
                beta_count = stats['beta_count']
                confidence_score = stats['average_confidence_score']
                chain_distance = stats['chain_dist_secondary_struct']
            
            # Write gene name and site information
                out_file.write(f"{gene_name}\t{residue_number}\t{stability_score}\t{helix_count}\t{beta_count}\t{confidence_score}\t\t{chain_distance}\t")
            
            # Write 'X' or '' for each database based on phosphorylation_info
                for database in base_databases:
                    site_info_list = phosphorylation_info.get(gene_name, [])
                    if any((database, True) in site_info for site_info in site_info_list):
                        out_file.write("X\t")
                    else:
                        out_file.write("-\t")

                out_file.write("\n")
