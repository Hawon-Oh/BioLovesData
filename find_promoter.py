# file path
# need to edit the path depending on your environment
fasta_file_path = "Bio/Semisulcospira habei HB4 DNA, Sha008, whole genome shotgun sequence.fasta"


# lenient validation of a DNA sequence based on specific criteria
def is_it_valid_seq(seq, min_seq_len, max_seq_len, allowed_nucleotides, max_ambiguity):

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
        # ambiguity = "number of un-allowed nucleotides" / len(seq)
        ambiguity = len(seq.translate(str.maketrans("","",allowed_nucleotides))) / len(seq)
        
        # return False if ambiguity exceeds the max limit
        if max_ambiguity < ambiguity:
            print(f"the ambiguity of the sequence exceeds the maximum: {max_ambiguity*100}% < {ambiguity*100}%")
            return False

    # return True if all checks are passed, indicating the sequence is valid
    return True

# ceil: round up the float into int
def ceil(value):
    if value > 0:  # when positive value
        integer_part = int(value)
        # it's floating number
        if value - integer_part > 0:
            return integer_part + 1  # add 1 to make it ceiling
        # no decimal place, just return integer
        else:
            return integer_part
    # when negative value, use int( ) and return
    elif value < 0:
        return int(value)
    # when value 0, return 0
    else:
        return 0


# using fuzz matching, find all subsequence in sequence
# and return two arrays: their locations & similarities
def find_all_subseq_in_seq(seq, subseq, threshold_float):
    locations = []           # list to store the positions where the subsequence is found
    similarities = []        # list to store the similarity value (0.0 ~ 1.0)
    subseq_len = len(subseq) # length of the subsequence
    seq_len = len(seq)       # length of the sequence

    # calculate the maximum int number of mismatches allowed based on the threshold
    maximum_unmatch_num = subseq_len - ceil(subseq_len * threshold_float)

    i = 0   # iteration counter
    iters = seq_len - subseq_len # number of iterations
    while i <= iters:
        unmatch_count = 0    # initialize the count of unmatched nucleotides

        # compare each nucleotide in the subsequence with the corresponding one in the main sequence
        for j in range(subseq_len):
             # when nucleotide at the current position does not match
            if seq[i + j] != subseq[j]:
                unmatch_count += 1
                # stop comparing if mismatches exceed the allowed maximum
                if unmatch_count > maximum_unmatch_num:
                    break

        # when number of unmathced nucleotides less or equal to max
        if unmatch_count <= maximum_unmatch_num:
            # store the starting position where the subsequence was found
            locations.append(i)
            # calculate the similarity as the proportion of matching nucleotides and store it to array
            similarities.append((subseq_len - unmatch_count) / subseq_len)

        i += 1  # move to the next position in the main sequence

    # return the locations and similarities of all found subsequences
    return locations, similarities



#------------------------------------------------------------------------------------------
# MAIN FUNCTION
#------------------------------------------------------------------------------------------

# user choice: read file
user_mode_choice =  bool(int(input("Will you input sequence?(then input 0), or will read file?(then input 1)")))
if user_mode_choice:
    # open file in read mode
    file = open(fasta_file_path, 'r')  # fasta_file_path is definded in the line 1
    # get info
    fasta_description = file.readline()
    seq = file.read().replace('\n', '')
    # close file
    file.close()
    # print description
    print(fasta_description)

# user choice: input sequence
else:
    seq = str(input("Input the sequence: "))


# get user input: length of sequence
print(f"Sequence length is: {len(seq):,}")
user_seq_len = int(input("How long length you want to analyze?: "))
this_seq = seq[:user_seq_len] # 10,000,000

# get user input: maximum ambiguity
max_ambiguity = float(input("Maximum ambiguity? (0.001 ~ 1.0): "))

# check validation
allowed_nucleotides = "ATGC"
min_seq_length = 28          # promoter_len * 2 + min_spacer_len
max_seq_length = 500000000   # 500,000,000
is_it_valid_seq = is_it_valid_seq(this_seq, min_seq_length, max_seq_length, allowed_nucleotides, max_ambiguity)


# If the sequence is valid, proceed to find locations and similarities for -10 and -35 boxes
box_10 = "TATAAT"   # pribnow_box sequence
box_35 = "TTGACA"   # -35 box sequence
threshold = 0.6     # Similarity threshold (60% similarity)

if is_it_valid_seq:
    # Find all locations and similarities for the -10 box in the sequence
    box_10_locations, box_10_similarities = find_all_subseq_in_seq(this_seq, box_10, threshold)
    # Find all locations and similarities for the -35 box in the sequence
    box_35_locations, box_35_similarities = find_all_subseq_in_seq(this_seq, box_35, threshold)

    # Optimal Sequence Selection:
    # find all possible promoter regions
    min_spacer_len = 16               # minimum allowable distance between -10 and -35 boxes
    max_spacer_len = 19               # maximum allowable distance between -10 and -35 boxes
    promoter_region_locations = []    # to store promoter region locations
    promoter_region_similarities = [] # to store promoter region similarities

    # find all possible promoter regions
    i = 0     # Index for the current position in box_10_locations
    j = 0     #               ''               in box_35_locations

    # iterate through both lists of locations while both indices are within their respective bounds
    while i < len(box_10_locations) and j < len(box_35_locations):

        # if -35 box is further than possible maximum location -> get next box 10
        if box_10_locations[i] + len(box_10) + max_spacer_len < box_35_locations[j]:
            i+=1  # move to next box_10

        # if -35 box is closer than possible minimum location -> get next box 35
        elif box_10_locations[i] + len(box_10) + min_spacer_len > box_35_locations[j]:
            j+=1  # move to next box_35

        # when the -10 and -35 boxes are within the allowable spacer range
        else:
            # add the current pair of -10 and -35 box locations to the promoter region location array
            promoter_region_locations.append([box_10_locations[i], box_35_locations[j]])
            # add the corresponding similarities for the -10 and -35 boxes to the promoter region similarity array
            promoter_region_similarities.append([box_10_similarities[i], box_35_similarities[j]])
            i+=1  # move to next box_10
            j+=1  #      ''      box_35



    #----------------------------------------------------
    # OUTPUT FORMATTING
    #----------------------------------------------------

    # display the number of promoter regions found 
    # prompt the user to input how many to print
    print(f"\nNumber of promoter regions is: {len(promoter_region_locations):,}")
    user_promoter_region_num = int(input("How many promoters regions to print?: "))

     # ANSI color codes for terminal output
    blue = "\033[34m"
    white = "\033[0m"   # to reset it to default text color
    yellow = "\033[33m"

    # print the table header for promoter regions
    print(f"\n{white}     -10 box(.sim)  spacer     -35 box(.sim): promoter regions")
    print("-------------------------------------------------------------------------------")
    
    # loop through and print the user-specified number of promoter regions
    i = 0
    while i < user_promoter_region_num:
        # get the start positions of the -10 and -35 boxes
        box_10_begin = promoter_region_locations[i][0]
        box_35_begin = promoter_region_locations[i][1]
        
        # print promoter region informations:
        # print box_10 location, box_10 similarity, spacer length, box_35 location, box_35 similarity, and promoter region sequence
        # ex:   23(0.67)     19           48(0.83): TACGAT TTTTTATCTGCAAGGGGAT TTTACA
        print(f"{white}{box_10_begin:>12,}({promoter_region_similarities[i][0]:.2f})     {white}{box_35_begin - box_10_begin - 6} {white}{box_35_begin:>12,}({promoter_region_similarities[i][1]:.2f}): {blue}{seq[box_10_begin:box_10_begin+6]} {yellow}{seq[box_10_begin+6:box_35_begin]:<19} {blue}{seq[box_35_begin:box_35_begin+6]}")
        i+=1    # Move to the next promoter region


# If the sequence is not valid, print an error message
else:
   print("sequence not valid")

# reset terminal text color back to default
print(f"{white}")