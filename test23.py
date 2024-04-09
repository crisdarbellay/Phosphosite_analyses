def align_and_find_matching_subsequences(sequence1, sequence2):
    matching_subsequences = []
    
    # Parcourir la séquence 1
    i=0
    for j in range(i + 5, len(sequence1) + 1):
        subsequence = sequence1[i:j]
            
            # Rechercher la sous-séquence dans la séquence 2
        if subsequence in sequence2:
            start_index = sequence2.index(subsequence)
            end_index = start_index + len(subsequence)
            counting = True
        else:
            if counting:
                counting=False
                matching_subsequences.append((i + 1, j-1))
            i=j
    
    # Formater les sous-séquences correspondantes
    formatted_subsequences = "/".join([f"A{start}-{end}" for start, end in matching_subsequences])
    
    return formatted_subsequences

# Exemple d'utilisation
sequence_pdb1 = "ADGSJFASUDSJBRJSABFJABASJBMM"
sequence_pdb2 = "BBAGADGSJFASUDSJTRJSABFJABASJB"

matching_subsequences = align_and_find_matching_subsequences(sequence_pdb1, sequence_pdb2)
print("Sous-séquences correspondantes dans la séquence PDB1 :", matching_subsequences)
