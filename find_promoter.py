# Modified version of findPromoter.py for web application

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

def is_valid_seq(seq, min_seq_len, max_seq_len, allowed_nucleotides, max_ambiguity):
    """Validate DNA sequence based on specific criteria"""
    if len(seq) < min_seq_len:
        return False, "Sequence length lesser than minimum length"
    elif len(seq) > max_seq_len:
        return False, "Sequence length longer than maximum length"
    else:
        ambiguity = len(seq.translate(str.maketrans("","",allowed_nucleotides))) / len(seq)
        if max_ambiguity < ambiguity:
            return False, f"The ambiguity of the sequence exceeds the maximum: {max_ambiguity*100}% < {ambiguity*100}%"
    return True, "Sequence is valid"

def ceil(value):
    """Custom ceiling function"""
    if value > 0:
        integer_part = int(value)
        if value - integer_part > 0:
            return integer_part + 1
        else:
            return integer_part
    elif value < 0:
        return int(value)
    else:
        return 0

def find_all_subseq_in_seq(seq, subseq, threshold_float):
    """Find all subsequences in sequence using fuzzy matching"""
    locations = []
    similarities = []
    subseq_len = len(subseq)
    seq_len = len(seq)
    
    maximum_unmatch_num = subseq_len - ceil(subseq_len * threshold_float)
    
    i = 0
    iters = seq_len - subseq_len
    while i <= iters:
        unmatch_count = 0
        
        for j in range(subseq_len):
            if seq[i + j] != subseq[j]:
                unmatch_count += 1
                if unmatch_count > maximum_unmatch_num:
                    break
        
        if unmatch_count <= maximum_unmatch_num:
            locations.append(i)
            similarities.append((subseq_len - unmatch_count) / subseq_len)
        
        i += 1
    
    return locations, similarities

def find_promoter_analysis(file_content, seq_length=None, max_ambiguity=0.001, max_results=10):
    """
    Main function for promoter analysis
    Args:
        file_content: String content of the uploaded file
        seq_length: Length of sequence to analyze (optional)
        max_ambiguity: Maximum ambiguity allowed (default: 0.001)
        max_results: Maximum number of results to return (default: 10)
    Returns:
        dict: Analysis results or error message
    """
    try:
        # Try to read as FASTA file first
        info_arr, seq_arr = read_seq_file_content(file_content)
        
        if seq_arr:
            # Use first sequence from FASTA file
            seq = seq_arr[0]
            description = info_arr[0] if info_arr else "Sequence from uploaded file"
        else:
            # Treat as raw sequence
            seq = file_content.replace('\n', '').replace(' ', '').upper()
            description = "Raw sequence input"
        
        if not seq:
            return {"error": "No sequence found in the input"}
        
        # Use provided length or full sequence
        if seq_length and seq_length < len(seq):
            this_seq = seq[:seq_length]
        else:
            this_seq = seq
        
        # Validation parameters
        allowed_nucleotides = "ATGC"
        min_seq_length = 28
        max_seq_length = 500000000
        
        is_valid, message = is_valid_seq(this_seq, min_seq_length, max_seq_length, allowed_nucleotides, max_ambiguity)
        
        if not is_valid:
            return {"error": message}
        
        # Promoter search parameters
        box_10 = "TATAAT"
        box_35 = "TTGACA"
        threshold = 0.6
        min_spacer_len = 16
        max_spacer_len = 19
        
        # Find -10 and -35 boxes
        box_10_locations, box_10_similarities = find_all_subseq_in_seq(this_seq, box_10, threshold)
        box_35_locations, box_35_similarities = find_all_subseq_in_seq(this_seq, box_35, threshold)
        
        # Find promoter regions
        promoter_regions = []
        
        i = 0
        j = 0
        while i < len(box_10_locations) and j < len(box_35_locations):
            if box_10_locations[i] + len(box_10) + max_spacer_len < box_35_locations[j]:
                i += 1
            elif box_10_locations[i] + len(box_10) + min_spacer_len > box_35_locations[j]:
                j += 1
            else:
                box_10_begin = box_10_locations[i]
                box_35_begin = box_35_locations[j]
                spacer_length = box_35_begin - box_10_begin - 6
                
                promoter_region = {
                    "box_10_location": box_10_begin,
                    "box_10_similarity": round(box_10_similarities[i], 2),
                    "box_35_location": box_35_begin,
                    "box_35_similarity": round(box_35_similarities[j], 2),
                    "spacer_length": spacer_length,
                    "box_10_sequence": this_seq[box_10_begin:box_10_begin+6],
                    "spacer_sequence": this_seq[box_10_begin+6:box_35_begin],
                    "box_35_sequence": this_seq[box_35_begin:box_35_begin+6],
                    "full_promoter": this_seq[box_10_begin:box_35_begin+6]
                }
                promoter_regions.append(promoter_region)
                i += 1
                j += 1
        
        # Limit results
        if len(promoter_regions) > max_results:
            promoter_regions = promoter_regions[:max_results]
        
        return {
            "success": True,
            "description": description,
            "sequence_length": len(this_seq),
            "total_promoters": len(promoter_regions),
            "promoter_regions": promoter_regions,
            "parameters": {
                "box_10": box_10,
                "box_35": box_35,
                "threshold": threshold,
                "min_spacer_len": min_spacer_len,
                "max_spacer_len": max_spacer_len,
                "max_ambiguity": max_ambiguity
            }
        }
        
    except Exception as e:
        return {"error": f"Analysis failed: {str(e)}"}

