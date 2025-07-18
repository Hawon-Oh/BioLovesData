# Modified version of consensusSymbols.py for web application

# Grantham distance scores
GRANTHAM_DISTANCE_SCORES = {
    "A": {"A": 0, "R": 58, "N": 64, "D": 89, "C": 92, "Q": 112, "E": 107, "G": 60, "H": 148, "I": 94, "L": 96, "K": 106, "M": 84, "F": 113, "P": 97, "S": 148, "T": 69, "W": 148, "Y": 112, "V": 88},
    "R": {"A": 58, "R": 0, "N": 89, "D": 121, "C": 144, "Q": 112, "E": 65, "G": 124, "H": 112, "I": 56, "L": 110, "K": 46, "M": 80, "F": 142, "P": 145, "S": 177, "T": 74, "W": 135, "Y": 155, "V": 99},
    "N": {"A": 64, "R": 89, "N": 0, "D": 46, "C": 139, "Q": 53, "E": 42, "G": 94, "H": 68, "I": 153, "L": 149, "K": 101, "M": 142, "F": 158, "P": 91, "S": 46, "T": 65, "W": 174, "Y": 143, "V": 133},
    "D": {"A": 89, "R": 121, "N": 46, "D": 0, "C": 154, "Q": 61, "E": 45, "G": 85, "H": 81, "I": 168, "L": 154, "K": 56, "M": 160, "F": 172, "P": 105, "S": 65, "T": 85, "W": 181, "Y": 160, "V": 152},
    "C": {"A": 92, "R": 144, "N": 139, "D": 154, "C": 0, "Q": 154, "E": 170, "G": 159, "H": 174, "I": 198, "L": 198, "K": 173, "M": 196, "F": 205, "P": 169, "S": 112, "T": 149, "W": 215, "Y": 194, "V": 192},
    "Q": {"A": 112, "R": 112, "N": 53, "D": 61, "C": 154, "Q": 0, "E": 29, "G": 87, "H": 68, "I": 109, "L": 113, "K": 53, "M": 101, "F": 116, "P": 76, "S": 68, "T": 89, "W": 130, "Y": 99, "V": 96},
    "E": {"A": 107, "R": 65, "N": 42, "D": 45, "C": 170, "Q": 29, "E": 0, "G": 98, "H": 74, "I": 134, "L": 138, "K": 56, "M": 119, "F": 140, "P": 93, "S": 80, "T": 107, "W": 152, "Y": 122, "V": 96},
    "G": {"A": 60, "R": 124, "N": 94, "D": 85, "C": 159, "Q": 87, "E": 98, "G": 0, "H": 98, "I": 184, "L": 184, "K": 127, "M": 184, "F": 192, "P": 42, "S": 56, "T": 59, "W": 184, "Y": 184, "V": 138},
    "H": {"A": 148, "R": 112, "N": 68, "D": 81, "C": 174, "Q": 68, "E": 74, "G": 98, "H": 0, "I": 94, "L": 98, "K": 24, "M": 87, "F": 100, "P": 77, "S": 89, "T": 47, "W": 115, "Y": 83, "V": 84},
    "I": {"A": 94, "R": 56, "N": 153, "D": 168, "C": 198, "Q": 109, "E": 134, "G": 184, "H": 94, "I": 0, "L": 5, "K": 102, "M": 10, "F": 21, "P": 95, "S": 142, "T": 89, "W": 61, "Y": 40, "V": 29},
    "L": {"A": 96, "R": 110, "N": 149, "D": 154, "C": 198, "Q": 113, "E": 138, "G": 184, "H": 98, "I": 5, "L": 0, "K": 107, "M": 15, "F": 22, "P": 98, "S": 153, "T": 92, "W": 61, "Y": 40, "V": 32},
    "K": {"A": 106, "R": 46, "N": 101, "D": 56, "C": 173, "Q": 53, "E": 56, "G": 127, "H": 24, "I": 102, "L": 107, "K": 0, "M": 95, "F": 102, "P": 103, "S": 121, "T": 78, "W": 110, "Y": 85, "V": 97},
    "M": {"A": 84, "R": 80, "N": 142, "D": 160, "C": 196, "Q": 101, "E": 119, "G": 184, "H": 87, "I": 10, "L": 15, "K": 95, "M": 0, "F": 28, "P": 87, "S": 135, "T": 81, "W": 67, "Y": 36, "V": 21},
    "F": {"A": 113, "R": 142, "N": 158, "D": 172, "C": 205, "Q": 116, "E": 140, "G": 192, "H": 100, "I": 21, "L": 22, "K": 102, "M": 28, "F": 0, "P": 114, "S": 155, "T": 103, "W": 37, "Y": 22, "V": 50},
    "P": {"A": 97, "R": 145, "N": 91, "D": 105, "C": 169, "Q": 76, "E": 93, "G": 42, "H": 77, "I": 95, "L": 98, "K": 103, "M": 87, "F": 114, "P": 0, "S": 74, "T": 38, "W": 147, "Y": 99, "V": 68},
    "S": {"A": 148, "R": 177, "N": 46, "D": 65, "C": 112, "Q": 68, "E": 80, "G": 56, "H": 89, "I": 142, "L": 153, "K": 121, "M": 135, "F": 155, "P": 74, "S": 0, "T": 58, "W": 177, "Y": 144, "V": 124},
    "T": {"A": 69, "R": 74, "N": 65, "D": 85, "C": 149, "Q": 89, "E": 107, "G": 59, "H": 47, "I": 89, "L": 92, "K": 78, "M": 81, "F": 103, "P": 38, "S": 58, "T": 0, "W": 148, "Y": 92, "V": 69},
    "W": {"A": 148, "R": 135, "N": 174, "D": 181, "C": 215, "Q": 130, "E": 152, "G": 184, "H": 115, "I": 61, "L": 61, "K": 110, "M": 67, "F": 37, "P": 147, "S": 177, "T": 148, "W": 0, "Y": 37, "V": 88},
    "Y": {"A": 112, "R": 155, "N": 143, "D": 160, "C": 194, "Q": 99, "E": 122, "G": 184, "H": 83, "I": 40, "L": 40, "K": 85, "M": 36, "F": 22, "P": 99, "S": 144, "T": 92, "W": 37, "Y": 0, "V": 55},
    "V": {"A": 88, "R": 99, "N": 133, "D": 152, "C": 192, "Q": 96, "E": 96, "G": 138, "H": 84, "I": 29, "L": 32, "K": 97, "M": 21, "F": 50, "P": 68, "S": 124, "T": 69, "W": 88, "Y": 55, "V": 0},
}

def consensus_symbols_analysis(sequences_input):
    """
    Analyze consensus symbols for multiple sequence alignments
    Args:
        sequences_input: List of aligned sequences or string with space-separated sequences
    Returns:
        dict: Analysis results with sequences and symbols
    """
    try:
        # Handle input format
        if isinstance(sequences_input, str):
            seqs = sequences_input.split()
        elif isinstance(sequences_input, list):
            seqs = sequences_input
        else:
            return {"error": "Invalid input format. Expected string or list of sequences."}
        
        if not seqs:
            return {"error": "No sequences provided"}
        
        # Check if all sequences have the same length
        seq_lengths = [len(seq) for seq in seqs]
        if len(set(seq_lengths)) != 1:
            return {"error": "All sequences must have the same length for alignment analysis"}
        
        # Constant symbols
        VARIABLE = '.'
        CONSERVATIVE = ':'
        GAP = ' '
        IDENTICAL = '*'
        
        symbols = [IDENTICAL] * len(seqs[0])
        allowedAcids = GRANTHAM_DISTANCE_SCORES.keys()
        
        # Iterate through sites
        for i in range(len(seqs[0])):
            # Iterate from 0 to last with picking pair of 2 adjacent seqs
            for j in range(0, len(seqs)-1, 2):
                if j+1 >= len(seqs):
                    break
                    
                seq1 = seqs[j]
                seq2 = seqs[j+1]
                
                # If any site is a gap
                if seq1[i] not in allowedAcids or seq2[i] not in allowedAcids:
                    symbols[i] = GAP
                    break
                
                # Variable - when score is more or equal to 100
                if GRANTHAM_DISTANCE_SCORES[seq1[i]][seq2[i]] >= 100:
                    symbols[i] = VARIABLE
                
                # Conservative - when any site isn't variable
                elif GRANTHAM_DISTANCE_SCORES[seq1[i]][seq2[i]] > 0 and symbols[i] != VARIABLE:
                    symbols[i] = CONSERVATIVE
        
        # Count symbol statistics
        symbol_counts = {
            'identical': symbols.count(IDENTICAL),
            'conservative': symbols.count(CONSERVATIVE),
            'variable': symbols.count(VARIABLE),
            'gap': symbols.count(GAP)
        }
        
        return {
            "success": True,
            "sequences": seqs,
            "symbols": ''.join(symbols),
            "symbol_counts": symbol_counts,
            "total_sites": len(symbols),
            "legend": {
                "*": "Identical",
                ":": "Conservative",
                ".": "Variable",
                " ": "Gap"
            }
        }
        
    except Exception as e:
        return {"error": f"Analysis failed: {str(e)}"}

