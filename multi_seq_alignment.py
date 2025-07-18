# Modified version of multiSeqAlignment.py for web application

import numpy as np

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
    
    lines = file_content.strip().split("\n")
    for line in lines:
        line = line.strip()
        
        if line.startswith(">"):
            if seq:
                seq_arr.append(seq)
                seq = ""
            info_arr.append(line[1:])
        else:
            seq += line
    
    if seq:
        seq_arr.append(seq)
    
    return info_arr, seq_arr

# BLOSUM62 Matrix
BLOSUM62 = {
    'A': {'A': 4,  'R': -1, 'N': -2, 'D': -2, 'C':  0, 'Q': -1, 'E': -1, 'G':  0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S':  1, 'T':  0, 'W': -3, 'Y': -2, 'V':  0},
    'R': {'A': -1, 'R':  5, 'N':  0, 'D': -2, 'C': -3, 'Q':  1, 'E':  0, 'G': -2, 'H':  0, 'I': -3, 'L': -2, 'K':  2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
    'N': {'A': -2, 'R':  0, 'N':  6, 'D':  1, 'C': -3, 'Q':  0, 'E':  0, 'G':  0, 'H':  1, 'I': -3, 'L': -3, 'K':  0, 'M': -2, 'F': -3, 'P': -2, 'S':  1, 'T':  0, 'W': -4, 'Y': -2, 'V': -3},
    'D': {'A': -2, 'R': -2, 'N':  1, 'D':  6, 'C': -3, 'Q':  0, 'E':  2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
    'C': {'A':  0, 'R': -3, 'N': -3, 'D': -3, 'C':  9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
    'Q': {'A': -1, 'R':  1, 'N':  0, 'D':  0, 'C': -3, 'Q':  5, 'E':  2, 'G': -2, 'H':  0, 'I': -3, 'L': -2, 'K':  1, 'M':  0, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
    'E': {'A': -1, 'R':  0, 'N':  0, 'D':  2, 'C': -4, 'Q':  2, 'E':  5, 'G': -2, 'H':  0, 'I': -3, 'L': -3, 'K':  1, 'M':  0, 'F': -3, 'P': -1, 'S':  0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'G': {'A':  0, 'R': -2, 'N':  0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G':  6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S':  0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
    'H': {'A': -2, 'R':  0, 'N':  1, 'D': -1, 'C': -3, 'Q':  0, 'E':  0, 'G': -2, 'H':  8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y':  2, 'V': -3},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I':  4, 'L':  2, 'K': -3, 'M':  1, 'F':  0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V':  3},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I':  2, 'L':  4, 'K': -2, 'M':  2, 'F':  0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V':  1},
    'K': {'A': -1, 'R': 2,  'N': 0,  'D': -1, 'C': -3, 'Q': 1,  'E': 1,  'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5,  'M': -1, 'F': -3, 'P': -1, 'S': 0,  'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0,  'E': 0,  'G': -3, 'H': -2, 'I': 1,  'L': 2,  'K': -1, 'M': 5,  'F': 0,  'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0,  'L': 0,  'K': -3, 'M': 0,  'F': 6,  'P': -4, 'S': -2, 'T': -2, 'W': 1,  'Y': 3,  'V': -1},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7,  'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
    'S': {'A': 1,  'R': -1, 'N': 1,  'D': 0,  'C': -1, 'Q': 0,  'E': 0,  'G': 0,  'H': -1, 'I': -2, 'L': -2, 'K': 0,  'M': -1, 'F': -2, 'P': -1, 'S': 4,  'T': 1,  'W': -3, 'Y': -2, 'V': -2},
    'T': {'A': 0,  'R': -1, 'N': 0,  'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1,  'T': 5,  'W': -2, 'Y': -2, 'V': 0},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1,  'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2,  'V': -3},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2,  'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3,  'P': -3, 'S': -2, 'T': -2, 'W': 2,  'Y': 7,  'V': -1},
    'V': {'A': 0,  'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3,  'L': 1,  'K': -2, 'M': 1,  'F': -1, 'P': -2, 'S': -2, 'T': 0,  'W': -3, 'Y': -1, 'V': 4}
}

TEAM_SAMWISE = {
    'A': {'A':  3, 'C':  0, 'D': -2, 'E':  0, 'F':  0, 'G':  1, 'H': -1, 'I':  0, 'K': -1, 'L': -1, 'M': -2, 'N': -1, 'P':  1, 'Q': -1, 'R': -3, 'S':  0, 'T':  0, 'V':  0, 'W': -2, 'Y': -2},
    'C': {'A':  0, 'C':  4, 'D':  0, 'E':  0, 'F': -3, 'G': -1, 'H': -1, 'I':  1, 'K': -3, 'L':  0, 'M':  0, 'N': -1, 'P':  1, 'Q': -1, 'R': -1, 'S': -1, 'T': -2, 'V':  1, 'W':  1, 'Y': -1},
    'D': {'A': -2, 'C':  0, 'D':  4, 'E':  1, 'F': -4, 'G':  0, 'H':  0, 'I': -1, 'K': -1, 'L':  0, 'M': -3, 'N': -1, 'P': -1, 'Q':  0, 'R':  0, 'S':  0, 'T': -1, 'V':  0, 'W':  1, 'Y':  0},
    'E': {'A':  0, 'C':  0, 'D':  1, 'E':  3, 'F': -1, 'G':  0, 'H':  0, 'I':  0, 'K': -2, 'L':  0, 'M':  0, 'N': -1, 'P': -1, 'Q':  0, 'R': -1, 'S':  0, 'T': -2, 'V':  1, 'W': -3, 'Y':  0},
    'F': {'A':  0, 'C': -3, 'D': -4, 'E': -1, 'F':  4, 'G': -1, 'H': -2, 'I':  0, 'K':  1, 'L':  0, 'M': -1, 'N': -1, 'P': -1, 'Q': -2, 'R': -1, 'S':  0, 'T':  1, 'V': -2, 'W':  0, 'Y':  2},
    'G': {'A':  1, 'C': -1, 'D':  0, 'E':  0, 'F': -1, 'G':  3, 'H': -2, 'I': -1, 'K': -1, 'L': -1, 'M':  2, 'N': -2, 'P': -2, 'Q':  1, 'R': -1, 'S':  1, 'T':  0, 'V': -1, 'W':  0, 'Y':  0},
    'H': {'A': -1, 'C': -1, 'D':  0, 'E':  0, 'F': -2, 'G': -2, 'H':  4, 'I': -1, 'K': -1, 'L': -1, 'M': -4, 'N':  2, 'P':  1, 'Q':  2, 'R':  0, 'S': -1, 'T':  0, 'V': -1, 'W': -2, 'Y':  1},
    'I': {'A':  0, 'C':  1, 'D': -1, 'E':  0, 'F':  0, 'G': -1, 'H': -1, 'I':  3, 'K':  1, 'L':  0, 'M':  0, 'N': -2, 'P': -1, 'Q': -1, 'R':  1, 'S': -1, 'T':  0, 'V':  1, 'W': -2, 'Y': -1},
    'K': {'A': -1, 'C': -3, 'D': -1, 'E': -2, 'F':  1, 'G': -1, 'H': -1, 'I':  1, 'K':  2, 'L':  0, 'M':  1, 'N':  0, 'P': -1, 'Q':  0, 'R':  1, 'S': -1, 'T':  1, 'V':  0, 'W': -1, 'Y':  0},
    'L': {'A': -1, 'C':  0, 'D':  0, 'E':  0, 'F':  0, 'G': -1, 'H': -1, 'I':  0, 'K':  0, 'L':  2, 'M': -1, 'N':  0, 'P': -1, 'Q': -1, 'R': -1, 'S':  0, 'T':  0, 'V':  0, 'W': -2, 'Y':  0},
    'M': {'A': -2, 'C':  0, 'D': -3, 'E':  0, 'F': -1, 'G':  2, 'H': -4, 'I':  0, 'K':  1, 'L': -1, 'M':  5, 'N': -1, 'P': -3, 'Q':  2, 'R':  1, 'S':  0, 'T': -3, 'V': -1, 'W':  3, 'Y': -4},
    'N': {'A': -1, 'C': -1, 'D': -1, 'E': -1, 'F': -1, 'G': -2, 'H':  2, 'I': -2, 'K':  0, 'L':  0, 'M': -1, 'N':  3, 'P':  0, 'Q':  1, 'R': -1, 'S': -1, 'T': -2, 'V':  0, 'W': -1, 'Y':  0},
    'P': {'A':  1, 'C':  1, 'D': -1, 'E': -1, 'F': -1, 'G': -2, 'H':  1, 'I': -1, 'K': -1, 'L': -1, 'M': -3, 'N':  0, 'P':  4, 'Q':  0, 'R': -3, 'S': -1, 'T':  1, 'V':  0, 'W':  0, 'Y': -1},
    'Q': {'A': -1, 'C': -1, 'D':  0, 'E':  0, 'F': -2, 'G':  1, 'H':  2, 'I': -1, 'K':  0, 'L': -1, 'M':  2, 'N':  1, 'P':  0, 'Q':  3, 'R':  0, 'S': -1, 'T': -1, 'V':  0, 'W':  0, 'Y': -3},
    'R': {'A': -3, 'C': -1, 'D':  0, 'E': -1, 'F': -1, 'G': -1, 'H':  0, 'I':  1, 'K':  1, 'L': -1, 'M':  1, 'N': -1, 'P': -3, 'Q':  0, 'R':  4, 'S': -1, 'T':  0, 'V':  0, 'W':  0, 'Y': -1},
    'S': {'A':  0, 'C': -1, 'D':  0, 'E':  0, 'F':  0, 'G':  1, 'H': -1, 'I': -1, 'K': -1, 'L':  0, 'M':  0, 'N': -1, 'P': -1, 'Q': -1, 'R': -1, 'S':  3, 'T':  0, 'V':  0, 'W':  1, 'Y': -1},
    'T': {'A':  0, 'C': -2, 'D': -1, 'E': -2, 'F':  1, 'G':  0, 'H':  0, 'I':  0, 'K':  1, 'L':  0, 'M': -3, 'N': -2, 'P':  1, 'Q': -1, 'R':  0, 'S':  0, 'T':  3, 'V': -2, 'W': -1, 'Y':  0},
    'V': {'A':  0, 'C':  1, 'D':  0, 'E':  1, 'F': -2, 'G': -1, 'H': -1, 'I':  1, 'K':  0, 'L':  0, 'M': -1, 'N':  0, 'P':  0, 'Q':  0, 'R':  0, 'S':  0, 'T': -2, 'V':  2, 'W':  0, 'Y': -1},
    'W': {'A': -2, 'C':  1, 'D':  1, 'E': -3, 'F':  0, 'G':  0, 'H': -2, 'I': -2, 'K': -1, 'L': -2, 'M':  3, 'N': -1, 'P':  0, 'Q':  0, 'R':  0, 'S':  1, 'T': -1, 'V':  0, 'W':  7, 'Y': -2},
    'Y': {'A': -2, 'C': -1, 'D':  0, 'E':  0, 'F':  2, 'G':  0, 'H':  1, 'I': -1, 'K':  0, 'L':  0, 'M': -4, 'N':  0, 'P': -1, 'Q': -3, 'R': -1, 'S': -1, 'T':  0, 'V': -1, 'W': -2, 'Y': -4},
}

defaultScoringMatrix = {
    'A': {'A': 10, 'C': -1, 'G': -1, 'T': -1},
    'C': {'A': -1, 'C': 10, 'G': -1, 'T': -1},
    'G': {'A': -1, 'C': -1, 'G': 10, 'T': -1},
    'T': {'A': -1, 'C': -1, 'G': -1, 'T': 10}
}

HoxD55 = {
    'A': {'A':   91, 'C':  -90, 'G':  -25, 'T': -100},
    'C': {'A':  -90, 'C':  100, 'G': -100, 'T':  -25},
    'G': {'A':  -25, 'C': -100, 'G':  100, 'T':  -90},
    'T': {'A': -100, 'C':  -25, 'G':  -90, 'T':   91}
}

HoxD70 = {
    'A': {'A': 91, 'C': -114, 'G': -31, 'T': -123},
    'C': {'A': -114, 'C': 100, 'G': -125, 'T': -31},
    'G': {'A': -31, 'C': -125, 'G': 100, 'T': -114},
    'T': {'A': -123, 'C': -31, 'G': -114, 'T': 91}
}

SCORING_MATRICES = {
    "BLOSUM62": BLOSUM62,
    "TEAM_SAMWISE": TEAM_SAMWISE,
    "defaultScoringMatrix": defaultScoringMatrix,
    "HoxD55": HoxD55,
    "HoxD70": HoxD70
}

def seqAlign(seqLeft, seqUp, scoringMatrix, gap_penalty=-2):
    """Single sequence alignment function"""
    DIAGONAL = 1
    UP = 2
    LEFT = 3

    scoreTable = np.zeros((len(seqLeft)+1, len(seqUp)+1), dtype=np.int32)
    traceTable = scoreTable.copy()
    scoreTable[0] = list(range(0, (len(seqUp)+1) * gap_penalty, gap_penalty))
    scoreTable[:, 0] = list(range(0, (len(seqLeft)+1) * gap_penalty, gap_penalty))

    allowedAcids = scoringMatrix.keys()

    for i in range(1,len(seqLeft)+1):
        for j in range(1,len(seqUp)+1):
            MatchOrMismatch = max((scoringMatrix[_][seqUp[j-1]] for _ in seqLeft[i-1] if _ in allowedAcids and seqUp[j-1] in allowedAcids),default=gap_penalty)

            diagonalScore = scoreTable[i-1, j-1] + MatchOrMismatch
            upScore = scoreTable[i-1, j] + gap_penalty
            leftScore = scoreTable[i, j-1] + gap_penalty
            maxScore = max(diagonalScore, upScore, leftScore)
            scoreTable[i, j] = maxScore
            
            if maxScore == diagonalScore: 
                traceTable[i, j] = DIAGONAL
            elif maxScore == upScore:     
                traceTable[i, j] = UP
            else:                         
                traceTable[i, j] = LEFT

    alignedSeq = []
    mergedSeqs = []

    i = len(seqLeft)
    j = len(seqUp)
    while i > 0 or j > 0:
        if traceTable[i, j] == DIAGONAL:
            mergedSeqs.insert(0, ''.join(set(seqLeft[i-1] + seqUp[j-1])))
            alignedSeq.insert(0, seqUp[j-1])
            i -= 1
            j -= 1
        elif traceTable[i, j] == UP:
            mergedSeqs.insert(0, ''.join(set(seqLeft[i-1] + '-')))
            alignedSeq.insert(0, '-')
            i -= 1
        else:
            mergedSeqs.insert(0, ''.join(set(seqUp[j-1] + '-')))
            alignedSeq.insert(0, seqUp[j-1])
            j -= 1

    return alignedSeq, scoreTable[-1,-1], mergedSeqs

def multiSeqAlign(seqs, matrix_name, gap_penalty=-2):
    """Perform multiple sequence alignment"""
    if len(seqs) < 2:
        return None
    
    scoring_matrix = SCORING_MATRICES.get(matrix_name, BLOSUM62) # Default to BLOSUM62

    alignedSeqs = seqAlign(seqs[0], seqs[1], scoring_matrix, gap_penalty)[2]

    if len(seqs) > 2: 
        for seq in seqs[2:]:
            alignedSeqs = seqAlign(seqLeft=alignedSeqs, seqUp=seq, scoringMatrix=scoring_matrix, gap_penalty=gap_penalty)[2]

    alginedSeqTable = []
    for seq in seqs:
        alginedSeq = seqAlign(seqLeft=alignedSeqs, seqUp=seq, scoringMatrix=scoring_matrix, gap_penalty=gap_penalty)[0]
        alginedSeqTable.append(''.join(alginedSeq))

    return alginedSeqTable

def multi_seq_alignment_analysis(file_content, gap_penalty=-2, scoring_matrix_name="BLOSUM62"):
    """
    Main function for multiple sequence alignment analysis
    Args:
        file_content: String content of the uploaded file
        gap_penalty: Gap penalty value (default: -2)
        scoring_matrix_name: Name of the scoring matrix to use (default: "BLOSUM62")
    Returns:
        dict: Analysis results or error message
    """
    try:
        info_arr, seq_arr = read_seq_file_content(file_content)
        
        if not seq_arr:
            return {"error": "No sequences found in the file"}
        
        if len(seq_arr) < 2:
            return {"error": "At least 2 sequences are required for multiple sequence alignment"}
        
        # Perform multiple sequence alignment
        aligned_sequences = multiSeqAlign(seqs=seq_arr, matrix_name=scoring_matrix_name, gap_penalty=gap_penalty)
        
        if not aligned_sequences:
            return {"error": "Alignment failed"}
        
        # Create sequence information with names
        sequence_info = []
        for i, (name, original_seq, aligned_seq) in enumerate(zip(info_arr, seq_arr, aligned_sequences)):
            sequence_info.append({
                "id": i + 1,
                "name": name,
                "original_sequence": original_seq,
                "aligned_sequence": aligned_seq,
                "original_length": len(original_seq),
                "aligned_length": len(aligned_seq)
            })
        
        # Calculate alignment statistics
        alignment_length = len(aligned_sequences[0]) if aligned_sequences else 0
        gap_positions = []
        
        for pos in range(alignment_length):
            gap_count = sum(1 for seq in aligned_sequences if seq[pos] == '-')
            if gap_count > 0:
                gap_positions.append({
                    "position": pos + 1,
                    "gap_count": gap_count,
                    "gap_percentage": round((gap_count / len(aligned_sequences)) * 100, 2)
                })
        
        return {
            "success": True,
            "sequence_count": len(seq_arr),
            "alignment_length": alignment_length,
            "sequences": sequence_info,
            "aligned_sequences": aligned_sequences,
            "gap_statistics": {
                "total_gap_positions": len(gap_positions),
                "gap_positions": gap_positions[:20]  # Limit to first 20 for display
            },
            "parameters": {
                "scoring_matrix": scoring_matrix_name,
                "gap_penalty": gap_penalty
            }
        }
        
    except Exception as e:
        return {"error": f"Analysis failed: {str(e)}"}


