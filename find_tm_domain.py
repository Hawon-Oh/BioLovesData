# Modified version of findTmDomain.py for web application

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

def validate_seq(seq, min_seq_len, max_seq_len, allowed_nucleotides, max_ambiguity=0.0):
    """Validate sequence based on specific criteria"""
    if len(seq) < min_seq_len:
        return False, "Sequence length lesser than minimum length"
    elif len(seq) > max_seq_len:
        return False, "Sequence length longer than maximum length"
    else:
        ambiguity = len(seq.translate(str.maketrans("","",allowed_nucleotides))) / len(seq)
        if max_ambiguity < ambiguity:
            return False, f"The ambiguity of the sequence exceeds the maximum: {max_ambiguity*100}% < {ambiguity*100:.2f}%"
    return True, "Sequence is valid"

def get_scale_values(scale="KyteDoolittle"):
    """Return dictionary of hydrophobicity values for different scales"""
    
    # Kyte Doolittle
    kd_hydrophobicity = {
        'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
        'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
        'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
        'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
    }

    # Wimley White
    ww_hydrophobicity = {
        'I': 0.31, 'V': -0.07, 'L': 0.56, 'F': 1.13, 'C': 0.24,
        'M': 0.23, 'A': -0.17, 'G': -0.01, 'T': -0.14, 'S': -0.13,
        'W': 1.85, 'Y': 0.94, 'P': -0.45, 'H': -0.96, 'E': -2.02,
        'Q': -0.58, 'D': -1.23, 'N': -0.42, 'K': -0.99, 'R': -0.81
    }

    # Hessa Higy
    hh_hydrophobicity = {
        'I': -0.60, 'V': -0.31, 'L': -0.55, 'F': -0.32, 'C': -0.13,
        'M': -0.10, 'A': 0.11, 'G': 0.74, 'T': 0.52, 'S': 0.84,
        'W': 0.30, 'Y': 0.68, 'P': 2.23, 'H': 2.06, 'E': 2.68,
        'Q': 2.36, 'D': 3.49, 'N': 2.05, 'K': 2.71, 'R': 2.58
    }

    # Moon Fleming
    mf_hydrophobicity = {
        'I': -1.56, 'V': -0.78, 'L': -1.81, 'F': -2.20, 'C': 0.49,
        'M': -0.76, 'A': 0.00, 'G': 1.72, 'T': 1.78, 'S': 1.83,
        'W': -0.38, 'Y': -1.09, 'P': -1.52, 'H': 4.76, 'E': 1.64,
        'Q': 3.01, 'D': 2.95, 'N': 3.47, 'K': 5.39, 'R': 3.71
    }

    # Tanford Taylor
    tt_hydrophobicity = {
        'I': 1.97, 'V': 1.46, 'L': 1.82, 'F': 1.98, 'C': -0.30,
        'M': 1.40, 'A': 0.38, 'G': -0.19, 'T': -0.32, 'S': -0.53,
        'W': 1.53, 'Y': 0.49, 'P': -1.44, 'H': -1.44, 'E': -2.90,
        'Q': -1.84, 'D': -3.27, 'N': -1.62, 'K': -3.46, 'R': -2.57
    }

    hydrophobicity_scale_names = {
        'KyteDoolittle': kd_hydrophobicity, 'kd': kd_hydrophobicity,
        'WimleyWhite': ww_hydrophobicity,   'ww': ww_hydrophobicity,
        'HessaHigy': hh_hydrophobicity,     'hh': hh_hydrophobicity,
        'MoonFleming': mf_hydrophobicity,   'mf': mf_hydrophobicity,
        'TanfordTaylor': tt_hydrophobicity, 'tt': tt_hydrophobicity
    }

    if scale not in hydrophobicity_scale_names:
        raise ValueError("Invalid hydrophobicity scale.")

    return hydrophobicity_scale_names[scale]

def hydrophobicity(seq, scale="KyteDoolittle"):
    """Calculate the hydrophobicity of a sequence"""
    hydrophobicity_value = sum(get_scale_values(scale).get(aa, 0) for aa in seq) / len(seq)
    return hydrophobicity_value

def find_tm_domain(seq, scale="KyteDoolittle", min_window=18, max_window=21, threshold=1.5):
    """Find transmembrane domains in a sequence"""
    if len(seq) < min_window:
        raise ValueError("Sequence length is too short")
        
    domain_table = []

    for i in range(len(seq) - min_window + 1):
        for current_window in range(min_window, max_window + 1):
            hydroph_value = hydrophobicity(seq[i:i + current_window], scale)
            if hydroph_value >= threshold:
                domain_table.append({
                    "location": i,
                    "length": current_window,
                    "hydrophobicity": round(hydroph_value, 3),
                    "sequence": seq[i:i + current_window]
                })

    # Remove overlapping domains, keeping the one with higher hydrophobicity
    if len(domain_table) > 1:
        i = 0
        while i < len(domain_table) - 1:
            current = domain_table[i]
            next_domain = domain_table[i + 1]
    
            if current["location"] + current["length"] > next_domain["location"]:
                if current["hydrophobicity"] < next_domain["hydrophobicity"]:
                    del domain_table[i]
                else:
                    del domain_table[i + 1]
            else:
                i += 1

    return domain_table

def find_tm_domain_analysis(file_content, threshold=1.5, scale="KyteDoolittle"):
    """
    Main function for TM domain analysis
    Args:
        file_content: String content of the uploaded file
        threshold: Hydrophobicity threshold (default: 1.5)
        scale: Hydrophobicity scale to use (default: KyteDoolittle)
    Returns:
        dict: Analysis results or error message
    """
    try:
        # Try to read as FASTA file first
        info_arr, seq_arr = read_seq_file_content(file_content)
        
        if not seq_arr:
            # Treat as raw sequence
            seq_arr = [file_content.replace('\n', '').replace(' ', '').upper()]
            info_arr = ["Raw sequence input"]
        
        # Validation parameters
        residue_types = "IVLFCMAGTSWYPHEQDNKR"
        min_seq_len = 20
        max_seq_len = 500000000
        max_ambiguity = 0.02
        
        # TM domain parameters
        tm_domain_min_len = 18
        tm_domain_max_len = 21
        
        results = []
        
        for one_info, one_seq in zip(info_arr, seq_arr):
            is_valid, message = validate_seq(one_seq, min_seq_len, max_seq_len, residue_types, max_ambiguity)
            
            if not is_valid:
                results.append({
                    "name": one_info,
                    "sequence": one_seq,
                    "valid": False,
                    "error": message,
                    "tm_domains": []
                })
            else:
                tm_domain_table = find_tm_domain(
                    seq=one_seq, 
                    min_window=tm_domain_min_len, 
                    max_window=tm_domain_max_len, 
                    threshold=threshold, 
                    scale=scale
                )
                
                # Create sequence with TM domain annotations
                annotated_sequence = []
                last_pos = 0
                
                for domain in tm_domain_table:
                    # Add non-TM region
                    if domain["location"] > last_pos:
                        annotated_sequence.append({
                            "type": "normal",
                            "sequence": one_seq[last_pos:domain["location"]]
                        })
                    
                    # Add TM domain
                    annotated_sequence.append({
                        "type": "tm_domain",
                        "sequence": domain["sequence"],
                        "hydrophobicity": domain["hydrophobicity"]
                    })
                    
                    last_pos = domain["location"] + domain["length"]
                
                # Add remaining sequence
                if last_pos < len(one_seq):
                    annotated_sequence.append({
                        "type": "normal",
                        "sequence": one_seq[last_pos:]
                    })
                
                results.append({
                    "name": one_info,
                    "sequence": one_seq,
                    "valid": True,
                    "tm_domains": tm_domain_table,
                    "tm_domain_count": len(tm_domain_table),
                    "annotated_sequence": annotated_sequence
                })
        
        return {
            "success": True,
            "results": results,
            "parameters": {
                "threshold": threshold,
                "scale": scale,
                "min_window": tm_domain_min_len,
                "max_window": tm_domain_max_len
            },
            "available_scales": ["KyteDoolittle", "WimleyWhite", "HessaHigy", "MoonFleming", "TanfordTaylor"]
        }
        
    except Exception as e:
        return {"error": f"Analysis failed: {str(e)}"}

