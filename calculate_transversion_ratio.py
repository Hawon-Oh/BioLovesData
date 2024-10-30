
# THIS IS SAME FUNCTION IN CODE FOR PROBLEM 1
# read '.txt' or '.fasta' file
# return two arrays:
#   'info_arr': corresponding headers (names, descriptions)
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

                # store the info without '>'
                info_arr.append(line[1:])  
            
            # if not start with '>',
            # continue adding the sequence lines
            else:
                seq += line

        # add the last sequence after loop
        if seq:
            seq_arr.append(seq)

    # 'info_arr' contains all the sequence headers (names, descriptions)
    # 'seq_arr' contains the corresponding sequences
    return info_arr, seq_arr



# check if every sequence's length is same value
# validate_seq_len will raise ValueError if it's False. 
def validate_seq_len(seqs):
    len_arr = [len(seq) for seq in seqs]

    # if there is sequence length of different
    # set function erases duplicated value. ex. set([1, 1, 2, 2, 2]) = [1, 2]
    if len(set(len_arr)) != 1:
        print("One or more sequence has different length. Can't calculate transversion ratio.")
        return False
    # all sequence lengths are same. It's validate.
    return True



# check if two bases are a transversion
# this function is used in transversion_ratio function
def is_transversion(nucleotide1, nucleotide2):

    # if two nucleotides are same, not transversion
    if nucleotide1 == nucleotide2:
        return False

    purines = ['A', 'G']
    pyrimidines = ['C', 'T']
    
    # if first base = purine & second base= pyrimidine
    # or second base= purine & first base = pyrimidine,
    # it's transversion
    is_it_transversion = (nucleotide1 in pyrimidines and nucleotide2 in purines) or (nucleotide1 in purines and nucleotide2 in pyrimidines)
    # if transversion, return True
    return is_it_transversion



# check if two bases are a transition
# this function is used in transversion_ratio function
def is_transition(nucleotide1, nucleotide2):

    # if two nucleotides are same
    # no mutation, return false
    if nucleotide1 == nucleotide2:
        return False

    purines = ['A', 'G']
    pyrimidines = ['C', 'T']
    
    # if first base = purine & second base= pyrimidine
    # or second base= purine & first base = pyrimidine,
    # if transition, return True
    is_it_transition = (nucleotide1 in purines and nucleotide2 in purines) or (nucleotide1 in pyrimidines and nucleotide2 in pyrimidines)
    # if transition, return True
    return is_it_transition


# it's the way I used to solve (2-a).
# 1 - find most common nucleotide, comparing first nucleotides with most common one, count V & S
# 2 - loop to count V & S (compare a current nucleotide with next nucleotide)
def transversion_ratio(sequences):
    count_v = 0
    count_s = 0

    # find first elements from every sequence if not empty
    first_elements = [seq[0] for seq in sequences if seq]
    # calculate frequency and store it into dictionary
    frequency_dict = {ele: first_elements.count(ele) for ele in set(first_elements)}
    # find the most common nucleotides from first nucleotide of all sequences (I searched to learn MAX function.)
    most_common = max(frequency_dict, key=frequency_dict.get)

    # compare the most common nucleotide with each first element 
    # to count transversions, V and transitions, S
    for i in first_elements:
        if is_transversion(most_common, i): # if it's a transversion, count V
            count_v += 1
        elif is_transition(most_common, i): # if it's a transition, count S
            count_s += 1
        # if neither V or S, it means there was no mutation. so, don't count

    # loop through all sequences and compare consecutive nucleotides to count V and S
    for current_seq in range(len(sequences)):
        for current_nucleo in range(len(sequences[0]) - 1):
            nucleotide1 = sequences[current_seq][current_nucleo]       # current nucleotide
            nucleotide2 = sequences[current_seq][current_nucleo + 1]   # next nucleotide
            if is_transversion(nucleotide1, nucleotide2): # if it's a transversion, count V
                count_v += 1
            elif is_transition(nucleotide1, nucleotide2): # if it's a transition, count S
                count_s += 1                              # if not, count S
            # if neither V or S, it means there was no mutation. so, don't count

    # If transition(S) found, calculate & return the ratio
    if count_s != 0:
        ratio = count_v / count_s
    # If no transition(S) found, return infinity as the ratio
    else:
        ratio = float('inf')
    return count_v, count_s, ratio




#--------------------------------------------------------------
# MAIN FUNCTION
#--------------------------------------------------------------

filename = 'CS Midterm Problem 2.txt'
info_arr, seq_arr = read_seq_file(filename) # info_arr won't be used

if validate_seq_len(seq_arr):
    V, S, ratio = transversion_ratio(seq_arr)
    print(f'Transversions: {V}, Transitions: {S}, Ratio: {ratio}')
else:
    print("Sequence with invalid length found")
