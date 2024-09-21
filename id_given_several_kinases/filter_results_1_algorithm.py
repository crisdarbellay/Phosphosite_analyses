import os

def filter_lines_and_count(input_folder):
    subfolder_line_counts = {}
    
    # Iterate through each subfolder in the main folder
    for subfolder in os.listdir(input_folder):
        subfolder_path = os.path.join(input_folder, subfolder)
        if os.path.isdir(subfolder_path):
            total_lines_remaining = 0
            # Iterate through each file in the subfolder
            for file_name in os.listdir(subfolder_path):
                if file_name.endswith(".txt"):
                    input_file_path = os.path.join(subfolder_path, file_name)
                    
                    with open(input_file_path, 'r') as file:
                        lines = file.readlines()

                    header = lines[0]
                    filtered_lines = [header]
                    
                    for line in lines[1:]:
                        columns = line.strip().split('\t')
                        if len(columns) < 11:  # Ensure there are enough columns to check the conditions
                            continue
                        try:
                            col_10 = float(columns[10]) 
                            col_2 = float(columns[2]) 
                            col_11 = float(columns[11]) 
                        except:
                            continue
                        if col_10 < 70:
                            continue
                        if col_2 < 5 and col_11 > 1:
                            continue
                        if col_2 >= 5 or col_11 <= 1:
                            print(line.split()[0],line.split()[1])
                        
                        filtered_lines.append(line)

                    with open(input_file_path, 'w') as file:
                        file.writelines(filtered_lines)
                    
                    total_lines_remaining += len(filtered_lines) - 1  # Subtract 1 to exclude the header
            
            subfolder_line_counts[subfolder] = total_lines_remaining

    return subfolder_line_counts

# Example usage
input_folder = r"/mnt/c/Users/crisd/Desktop/PhosphoSwitch/datas/6_angstrom/control2/1_algorithm/all_kinases_results_filtered"
remaining_lines_counts = filter_lines_and_count(input_folder)
for subfolder, count in remaining_lines_counts.items():
    print(f"Total number of remaining lines in {subfolder}: {count}")