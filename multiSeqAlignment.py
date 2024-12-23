# seq files must be in the working directory, if no file, exit program
# this program will do multiple sequence alignment.

import os              # used for reading file
import numpy as np     # score&traceback matrix in seqAlign() are numpy array
import sys             # for exit() when no file found


#-------------------------------------------------------------------------
# READ FILE OR SEQUENCE FROM USER INPUT
#-------------------------------------------------------------------------

# read '.txt' or '.fasta' file
# return two arrays:
#   'name_arr': corresponding headers (names, descriptions)
#   'seq_arr' : one or more sequences
def read_seq_file(filename):
    seq_arr = []    # store one or more sequences
    info_arr = []   # store the corresponding headers (names/descriptions)
    seq = ""        # store one sequence that needs to be stored in 'seq_arr'
    
    # open the file in read mode
    with open(filename, 'r') as file:
        # read the file line by line
        for line in file:
            line = line.strip()  # remove trailing whitespace from the line
            
            # if line is the start of new sequence
            if line.startswith('>'):
                # if there's a previous sequence, store it
                if seq:
                    seq_arr.append(seq)
                    seq = ""     # reset sequence for next

                # store the name without '>'
                info_arr.append(line[1:])
                # name is first word in the header
                # name_arr.append(line[1:].split()[0]) header in test cases has space in each name, so can't use it

            # if not start with '>',
            # continue adding the sequence lines
            else:
                seq += line

        # add the last sequence after loop
        if seq:
            seq_arr.append(seq)

    # 'name_arr' contains all the sequence header
    # 'seq_arr' contains the corresponding (sequences)
    return info_arr, seq_arr


# print current directory path
print("Current directory:", os.getcwd())

# get text & fasta files in the current directory
seq_files = [f for f in os.listdir() if f.endswith('.txt') or f.endswith('.fasta')]

# if no file found, exit program
if not seq_files:
    print("No text or fasta files found in the current directory.")
    sys.exit(0) # exit

# if file found, print all files, so that user can choose a file
print("Readable files in the current directory:")
for order, file in enumerate(seq_files, start=1):
    print(f"{order} - {file}")

# iterate until user input is valid
# escape conditions are:
# user input is valid file number, valid filename, or sequence as string
while True:
    user_input = input("\nEnter a file number or filename: ")

    # if user input is file number, use the corresponding file 
    if user_input.isdigit():
        user_input = int(user_input) - 1 # -1 because array starts from 0
        # if valid file number, read&store file, and then break loop
        if 0 <= user_input < len(seq_files):
            names, seqs = read_seq_file(seq_files[user_input])
            break  # break loop.
        else:
            print("Invalid file number. Please try again.")

    # if user input is filename
    elif '.' in user_input:
        # if valid filename, read&store file, and then break loop
        if user_input in seq_files:
            names, seqs = read_seq_file(user_input)
            break # break loop
            # if unvalid filename, loop again
        else:
            print("Invalid file name. Please try again.")


# make seq dictionary with name (keys) $ seqs (values)
# dict in old python ver(<3.7) doesn't certain the order of values
# in this code, order doesn't matter, so old python is okay, just output order'll differ
seqTable = dict(zip(names, seqs))
#-------------------------------------------------------------------------





#-------------------------------------------------------------------------
# Scoring Matrices
#-------------------------------------------------------------------------

# default BLOSUM62 Matrix (from lecture slide)
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
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# Sequence Alignment
#-------------------------------------------------------------------------

def seqAlign(seqLeft, seqUp, scoringMatrix=BLOSUM62):
    """
    single sequence alignment function
    
    Args:
        seqLeft      : it can be normal seq (ex: ['M', 'T', 'S']) or merged seq(ex: ['M', '-T', '-SA', '-P', '-DQE])
                     : it's used to check similarity with seqUp
        seqUp        : normal seq (ex: ['M', 'T', 'S'])
                     : it'll be aligned based on seqLeft and returned
        scoringMatrix: dictionary table for scoring match/mismatch

    Returns:
        alignedSeq (char list): it's single aligned seq(param 'seqUp')
        scoreTable[-1,-1]     : similarity score between seqLeft & seqUp
        mergedSeqs            : (ex: ['M', '-T', '-SA', '-P', '-DQE]
                                it's merged seq from multiple seqs into one array.
                                it's used for multiple sequence alignment as param 'seqLeft'
    """
    # gap penalty. match score & mismatch penalty are included in scoringMatrix
    GAP_PENALTY = -2
    # represent directions in traceback
    DIAGONAL = 1
    UP = 2
    LEFT = 3

    # scoreTable & traceTable are used for calculating similarity score and traceback
    scoreTable = np.zeros((len(seqLeft)+1, len(seqUp)+1), dtype=np.int32)   # for storing all scores
    traceTable = scoreTable.copy()                                          # for traceback
    scoreTable[0] = list(range(0, (len(seqUp)+1) * GAP_PENALTY, GAP_PENALTY))               # initialize scores in 1st row
    scoreTable[:, 0] = list(range(0, (len(seqLeft)+1) * GAP_PENALTY, GAP_PENALTY))          #                   in 1st col

    # allowed amino acids array: it's used to exclude when wrong acid is given from sequence
    allowedAcids = scoringMatrix.keys()

    # iterate for every coordinate([i,j]) in Table (except 1st row and 1st col which initial gap penalties are applied)
    # and calculate score for each iteration
    for i in range(1,len(seqLeft)+1):
        for j in range(1,len(seqUp)+1):
            # first, calculate match/mismatch/gap penalty(or score)
            # this value will be added to score from diagonal direction

            # seqLeft can have multiple amino acids in one index
            # if seqLeft has multiple(merged) char(acid or gap) in current index, iterate for each char
            # and find the max score among them

            # if any of seqLeft or seqUp has unallowed acid, apply GAP penalty
            MatchOrMismatch = max((scoringMatrix[_][seqUp[j-1]] for _ in seqLeft[i-1] if _ in allowedAcids and seqUp[j-1] in allowedAcids),default=GAP_PENALTY)

            # calculate scores from 3 directions
            diagonalScore = scoreTable[i-1, j-1] + MatchOrMismatch
            upScore = scoreTable[i-1, j] + GAP_PENALTY
            leftScore = scoreTable[i, j-1] + GAP_PENALTY
            # store the highest score among 3 directions to scoreTable
            maxScore = max(diagonalScore, upScore, leftScore)
            scoreTable[i, j] = maxScore
            # store from which direction the highest score is, to traceTable
            if maxScore == diagonalScore: traceTable[i, j] = DIAGONAL
            elif maxScore == upScore:     traceTable[i, j] = UP
            else:                         traceTable[i, j] = LEFT


    # using traceback(traceTable), get an aligned sequence
    alignedSeq = [] # store aligned seqUp
    mergedSeqs = [] # store merged(seqLeft+seqUp) seq

    # traceback starts from the last index in traceTable, so get lengths
    i = len(seqLeft)
    j = len(seqUp)
    while i > 0 or j > 0:
        # when match/mismatch
        if traceTable[i, j] == DIAGONAL:
            # when mistmatch, add that mistmatched acid
            # mergedSeqs & alignedSeq is built from back to front
            mergedSeqs.insert(0, ''.join(set(seqLeft[i-1] + seqUp[j-1])))
            alignedSeq.insert(0, seqUp[j-1])
            i -= 1 # move to diagonal(up & left)
            j -= 1
        # when seqUp has gap
        elif traceTable[i, j] == UP:
            mergedSeqs.insert(0, ''.join(set(seqLeft[i-1] + '-')))
            alignedSeq.insert(0, '-')
            i -= 1 # move up
        # when seqLeft has gap
        else: # Left
            mergedSeqs.insert(0, ''.join(set(seqUp[j-1] + '-')))
            alignedSeq.insert(0, seqUp[j-1])
            j -= 1 # move to left

    return alignedSeq, scoreTable[-1,-1], mergedSeqs
#-------------------------------------------------------------------------



#-------------------------------------------------------------------------
# Multiple Seq Alignment
#-------------------------------------------------------------------------

# dynamically do multi Sequence alignments
def multiSeqAlign(seqs, matrix=BLOSUM62):

    # if there is no seq or only single seq, return
    if len(seqs) < 2:
        return None
    
    # do single seq alignment for 1st & 2nd sequences
    # mergedSeqs example: ['M', '-T', '-SA', '-P', '-DQE]
    alignedSeqs = seqAlign(seqs[0], seqs[1], matrix)[2]

    # if seqTable has more than 2 sequences, 
    # do multiple seq alignment iteratively for remain sequences
    if len(seqs) > 2: 
        for seq in seqs[2:]:
            alignedSeqs = seqAlign(seqLeft=alignedSeqs, seqUp=seq, scoringMatrix=matrix)[2]

    # using alignedSeqs, align every sequence one by one (dynamic alignment),
    # and store aligned sequence
    alginedSeqTable = []
    for seq in seqs:
        alginedSeq = seqAlign(seqLeft=alignedSeqs, seqUp=seq, scoringMatrix=matrix)[0] # do single alignment
        alginedSeqTable.append(''.join(alginedSeq))              # make list to str and append that str to list

    return alginedSeqTable # return list containing all aligned sequences
#-------------------------------------------------------------------------


# print result (no format isn't yet included)
result = multiSeqAlign(seqs=list(seqTable.values()), matrix=BLOSUM62)
for i in result:
    print(i, end="\n\n")

