import os


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



# lenient validation of a DNA sequence based on specific criteria
# return True or False
def validate_seq(seq, min_seq_len, max_seq_len, allowed_nucleotides, max_ambiguity=0.0):

    # return false if seq length lesser than minimum length
    if len(seq) < min_seq_len:
        print("seq length lesser than minimum length")  # print error message
        return False

    # return false if seq length longer than maximum length
    elif len(seq) > max_seq_len:
        print("seq length longer than maximum length")  # print error message
        return False

    # if length checks are passed, proceed to check for ambiguity
    else:
        # calculate ambiguity as the ratio of unallowed nucleotides to the total length of the sequence
        # ambiguity = "number of un-allowed nucleotides" / "length of seq"
        ambiguity = len(seq.translate(str.maketrans("","",allowed_nucleotides))) / len(seq)

        # return False if ambiguity exceeds the max limit
        if max_ambiguity < ambiguity:
            # print error message
            print(f"the ambiguity of the sequence exceeds the maximum: {max_ambiguity*100}% < {ambiguity*100 : .2f}%")
            return False

    # return True if all checks are passed, indicating the sequence is valid
    return True



# default scale is "KyteDoolittle"
# using short name('kd', 'ww', etc) is also available
# return dictionary of hydrophobicity values to use
def get_scale_values(scale="KyteDoolittle"):

    # I - Isoleucine
    # V - Valine
    # L - Leucine
    # F - Phenylalanine
    # C - Cysteine
    # M - Methionine
    # A - Alanine
    # G - Glycine
    # T - Threonine
    # S - Serine
    # W - Tryptophan
    # Y - Tyrosine
    # P - Proline
    # H - Histidine
    # E - Glutamic acid
    # Q - Glutamine
    # D - Aspartic acid
    # N - Asparagine
    # K - Lysine
    # R - Arginine
    
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

    # dictionary that maps scale names to their corresponding hydrophobicity values
    hydrophobicity_scale_names = {
        'KyteDoolittle': kd_hydrophobicity, 'kd': kd_hydrophobicity,
        'WimleyWhite': ww_hydrophobicity,   'ww': ww_hydrophobicity,
        'HessaHigy': hh_hydrophobicity,     'hh': hh_hydrophobicity,
        'MoonFleming': mf_hydrophobicity,   'mf': mf_hydrophobicity,
        'TanfordTaylor': tt_hydrophobicity, 'tt': tt_hydrophobicity
    }

    # check if the provided scale is valid
    if scale not in hydrophobicity_scale_names:
        raise ValueError("Invalid hydrophobicity scale.")

    # return hydrophobicity values as dictionary
    return hydrophobicity_scale_names[scale]



# calculate the hydrophobicity
# return float number
def hydrophobicity(seq, scale="KyteDoolittle"):
    # if invalid residue, use default value, 0 
    hydrophobicity_value = sum(get_scale_values(scale).get(aa, 0) for aa in seq) / len(seq)

    return hydrophobicity_value



# default scale is Kyte Doolittle
# to use Wimley White, try any of these:  scale="WimleyWhite"  |   scale="ww"
def find_tm_domain(seq, scale="KyteDoolittle", min_window=18, max_window=21, threshold=1.5):

    # check if sequence length is shorter than minimum window size,
    # this 'if' statement won't be used, because we use this after validating each sequence
    if len(seq) < min_window:
        raise ValueError("Sequence length is too short")
        
    domain_table = [] # to store found tm domains
                      # keys= 'loc'(location) | 'len'(length) | 'hydrophobicity')

    # loop over the starting position of the sequence
    for i in range(len(seq) - min_window + 1):
        # loop over possible lengths of the tm domain (window size)
        for current_window in range(min_window, max_window + 1):
            hydroph_value = hydrophobicity(seq[i:i + current_window], scale)
            if hydroph_value >= threshold:
                domain_table.append({"loc":i, "len":current_window, "hydrophobicity":hydroph_value})

    # if there are two or more domains
    # check any overlap to delete
    if len(domain_table) > 1:
        i = 0
        while i < len(domain_table) - 1:
            current = domain_table[i]
            next = domain_table[i + 1] # next tm domain to compare
    
            # if two domains are overlapped.
            # compare each their hydrophobicity
            # and delete domain with lower one.
            if current["loc"] + current["len"] > next["loc"]:
                # If the next domain's hydrophobicity is greater than the current one's
                # delete current domain
                # if not, delete next domain
                if current["hydrophobicity"] < next["hydrophobicity"]:
                    del domain_table[i]
                else:
                    del domain_table[i + 1]
            # if no overlap domain, increase index
            else:
                i += 1

    return domain_table # return the dictionary list of tm domains



#-------------------------------------------------------------------------
# MAIN PART STARTS HERE
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# READ FILE OR SEQUENCE FROM USER INPUT
#-------------------------------------------------------------------------
print("Current directory:", os.getcwd()) # print current directory path

# get text & fasta files in the current directory
seq_files = [f for f in os.listdir() if f.endswith('.txt') or f.endswith('.fasta')]

# if no file found
if not seq_files:
    print("No text or fasta files found in the current directory.")
# if file found, print all files as list
else:
    print("Readable files in the current directory:")
    for order, file in enumerate(seq_files, start=1):
        print(f"{order} - {file}")


# iterate until user input is valid
# escape conditions:
# user input is valid file number, valid filename, or sequence as string
while True:
    user_input = input("\nEnter a file number, filename, or sequence: ")

    # if user input is file number, use the corresponding file 
    if user_input.isdigit():
        user_input = int(user_input) - 1 # -1 because array starts from 0
        # if valid file number, read&store file, and then break loop
        if 0 <= user_input < len(seq_files):
            info_arr, seq_arr = read_seq_file(seq_files[user_input])
            break  # break loop.
        else:
            print("Invalid file number. Please try again.")

    # if user input is filename
    elif '.' in user_input:
        # if valid filename, read&store file, and then break loop
        if user_input in seq_files:
            info_arr, seq_arr = read_seq_file(user_input)
            break # break loop

        # if unvalid filename, loop again
        else:
            print("Invalid file name. Please try again.")

    # if user input is sequence
    else:
        info_arr = ["Sequence from user input"] # read header
        seq_arr = [user_input]                  # read sequence
        break  # user inputs the sequence. break loop



#-------------------------------------------------------------------------
# FIND TM DOMAIN FOR EACH SEQUENCE & PRINT WITH OUTPUT FORMAT
#-------------------------------------------------------------------------


# variables to validate sequence
residue_types = "IVLFCMAGTSWYPHEQDNKR" # legal residue types
min_seq_len = 20                       # adequate min length
max_seq_len = 500000000                # 500,000,000
max_ambiguity = 0.02                   # more illegal residue types, bigger ambiguity value

# variables to find tm domain
tm_domain_min_len = 18                 # will be used for minimum window size
tm_domain_max_len = 21                 # will be used for maximum window size
user_threshold = float(input("Enter threshold: ")) # get threshold from user input
user_scale = (input("Enter scale to use('KyteDoolittle','WimleyWhite', 'HessaHigy', 'MoonFleming', 'TanfordTaylor'): "))

# ANSI color codes for output format
blue = "\033[34m"   # explanation color
yellow = "\033[33m" # tm domain color
white = "\033[0m"   # default text color (& non-tm domain color)

# loop through each sequence and its corresponding information
# sequence array & information array have equal length
for one_info, one_seq in zip(info_arr, seq_arr):
    print(f"\n{blue}Name:{white} {one_info}")
    # if not valid sequence, print it's not valid
    if not validate_seq(one_seq, min_seq_len, max_seq_len, residue_types, max_ambiguity):
        print(f"Not valid sequence: {one_seq}")
        
    # if valid sequence, find tm domains and print them all with distinguishable color
    else: 
        # call function to find tm domains in sequence based on defined parameters
        tm_domain_table = find_tm_domain(seq=one_seq, min_window=tm_domain_min_len, max_window=tm_domain_max_len, threshold=user_threshold, scale=user_scale)

        # if there's no tm domain, print all sequence in white color
        if not tm_domain_table:
            print(f"{blue}Sequence:{white} {one_seq}", end="\n\n")

        # if there's tm domain, print dm domains in yellow color
        else:
            print(f"{blue}Sequence:{white} ", end="")
            
            # loop through the TM domain table to print each domain and surrounding sequences
            i = 0 # index to track current position in sequence
            for j in tm_domain_table:
                # print subsequence that isn't tm domain
                print(f"{one_seq[i:j['loc']]}", end="")
                # print tm domain in yellow color
                print(f"{yellow}{one_seq[j['loc']:j['loc']+j['len']]}{white}", end="")
                i = j['loc'] + j['len']     # update index to the end of current tm domain

            # after loop, print all remaining part of sequence
            print(f"{one_seq[i:]}", end="\n\n")