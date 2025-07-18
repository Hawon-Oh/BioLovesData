# Modified version of calculateTransversionRatio.py for web application

def read_seq_file_content(file_content):
    """
    Read sequence file content and return sequences
    Args:
        file_content: String content of the file
    Returns:
        info_arr: List of sequence headers
        seq_arr: List of sequences
    """
    seq_arr = []
    info_arr = []
    seq = ""
    
    lines = file_content.strip().split('\n')
    for line in lines:
        line = line.strip()
        
        if line.startswith('>'):
            if seq:
                seq_arr.append(seq)
                seq = ""
            info_arr.append(line[1:])
        else:
            seq += line
    
    if seq:
        seq_arr.append(seq)
    
    return info_arr, seq_arr

def validate_seq_len(seqs):
    """Check if every sequence's length is same value"""
    len_arr = [len(seq) for seq in seqs]
    
    if len(set(len_arr)) != 1:
        return False, "One or more sequence has different length. Can't calculate transversion ratio."
    return True, "All sequences have same length."

def is_transversion(nucleotide1, nucleotide2):
    """Check if two bases are a transversion"""
    if nucleotide1 == nucleotide2:
        return False

    purines = ['A', 'G']
    pyrimidines = ['C', 'T']
    
    is_it_transversion = (nucleotide1 in pyrimidines and nucleotide2 in purines) or (nucleotide1 in purines and nucleotide2 in pyrimidines)
    return is_it_transversion

def is_transition(nucleotide1, nucleotide2):
    """Check if two bases are a transition"""
    if nucleotide1 == nucleotide2:
        return False

    purines = ['A', 'G']
    pyrimidines = ['C', 'T']
    
    is_it_transition = (nucleotide1 in purines and nucleotide2 in purines) or (nucleotide1 in pyrimidines and nucleotide2 in pyrimidines)
    return is_it_transition

def transversion_ratio(sequences):
    """Calculate transversion ratio"""
    count_v = 0
    count_s = 0

    # Find first elements from every sequence if not empty
    first_elements = [seq[0] for seq in sequences if seq]
    # Calculate frequency and store it into dictionary
    frequency_dict = {ele: first_elements.count(ele) for ele in set(first_elements)}
    # Find the most common nucleotides from first nucleotide of all sequences
    most_common = max(frequency_dict, key=frequency_dict.get)

    # Compare the most common nucleotide with each first element 
    for i in first_elements:
        if is_transversion(most_common, i):
            count_v += 1
        elif is_transition(most_common, i):
            count_s += 1

    # Loop through all sequences and compare consecutive nucleotides
    for current_seq in range(len(sequences)):
        for current_nucleo in range(len(sequences[0]) - 1):
            nucleotide1 = sequences[current_seq][current_nucleo]
            nucleotide2 = sequences[current_seq][current_nucleo + 1]
            if is_transversion(nucleotide1, nucleotide2):
                count_v += 1
            elif is_transition(nucleotide1, nucleotide2):
                count_s += 1

    # Calculate ratio
    if count_s != 0:
        ratio = count_v / count_s
    else:
        ratio = float('inf')
    
    return count_v, count_s, ratio

def calculate_transversion_ratio_analysis(file_content):
    """
    Main function for transversion ratio analysis
    Args:
        file_content: String content of the uploaded file
    Returns:
        dict: Analysis results or error message
    """
    try:
        info_arr, seq_arr = read_seq_file_content(file_content)
        
        if not seq_arr:
            return {"error": "No sequences found in the file"}
        
        is_valid, message = validate_seq_len(seq_arr)
        
        if is_valid:
            V, S, ratio = transversion_ratio(seq_arr)
            return {
                "success": True,
                "transversions": V,
                "transitions": S,
                "ratio": ratio if ratio != float('inf') else "Infinity",
                "sequences_count": len(seq_arr),
                "message": f"Transversions: {V}, Transitions: {S}, Ratio: {ratio if ratio != float('inf') else 'Infinity'}"
            }
        else:
            return {"error": message}
            
    except Exception as e:
        return {"error": f"Analysis failed: {str(e)}"}

