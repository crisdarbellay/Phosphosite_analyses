import pandas as pd
from prettytable import PrettyTable

def sort_key_with_nan(val):
    return float('-inf') if pd.isna(val) else val

def write_sorted_results(file_path, header, results_sorted, use_pretty_table=False):
    if use_pretty_table:
        table = PrettyTable()
        table.field_names = header
        for result in results_sorted:
            table.add_row(result)
        table.border = False
        table.header = True
        with open(file_path, 'w') as f_out_sorted:
            f_out_sorted.write(table.get_string())
    else:
        with open(file_path, 'w') as f_out_sorted:
            f_out_sorted.write('\t'.join(header) + '\n')
            for result in results_sorted:
                f_out_sorted.write('\t'.join(map(str, result)) + '\n')

def process_file(input_file, output_txt, output_csv, output_file_sorted_d, output_file_sorted_e, output_file_sorted_e_diff, output_file_sorted_d_diff):
    df = pd.read_csv(input_file, sep='\t')
    
    header = df.columns.tolist()
    results = df.values.tolist()
    
    # Sorting and writing results
    results_sorted_d = sorted(results, key=lambda x: sort_key_with_nan(x[-2]), reverse=True)
    write_sorted_results(output_file_sorted_d, header, results_sorted_d, use_pretty_table=True)

    results_sorted_e = sorted(results, key=lambda x: sort_key_with_nan(x[-1]), reverse=True)
    write_sorted_results(output_file_sorted_e, header, results_sorted_e, use_pretty_table=True)

    results_sorted_e_diff = sorted(results, key=lambda x: sort_key_with_nan(x[-11]), reverse=True)
    write_sorted_results(output_file_sorted_e_diff, header, results_sorted_e_diff, use_pretty_table=True)

    results_sorted_d_diff = sorted(results, key=lambda x: sort_key_with_nan(x[-12]), reverse=True)
    write_sorted_results(output_file_sorted_d_diff, header, results_sorted_d_diff, use_pretty_table=True)



# Example usage
input_file = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final_full_contigs/results_rmsd.txt"
temp_file = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/temp_results.txt"
output_file = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final_full_contigs/results_rmsd.txt"
output_total = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final_full_contigs/results_total.txt"
output_file_sorted_d = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final_full_contigs/results_rmsd_sorted_d.txt"
output_file_sorted_e = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final_full_contigs/results_rmsd_sorted_e.txt"
output_file_sorted_d_diff = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final_full_contigs/results_rmsd_sorted_wt-60.txt"
output_file_sorted_e_diff = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final_full_contigs/results_rmsd_sorted_wt-100.txt"
output_txt = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final_full_contigs/results_table.txt"
output_csv = "/mnt/c/Users/crisd/Desktop/ProteinDesign/results_new/human/results_final_full_contigs/results_table.csv"

process_file(input_file, output_txt, output_csv, output_file_sorted_d, output_file_sorted_e, output_file_sorted_e_diff, output_file_sorted_d_diff)
