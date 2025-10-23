# Libraries
import re

# connect to github? ====

## Read in data ================================================================================
# open file
with open('ref-for-bwt-search.fa', 'r') as file: 
    # read separate lines
    line = file.readlines()
    # select out sequence-specific data
    seq = line[1] 
    

## Could check sequence for suitibility here - basically that its a compatible genetic code    

## Preparing the sequence for BWT =======================================================================

# make test sequence
seq = 'AGTCTATCTGA$'

# append special character and overwrite
#seq = seq + '$'

## Generate BWT ==========================================================================
# set up empty list
bwt = []
# identify the last character
last_character = ''

# generate rotations and append to matrix
while last_character != '$': 
    # add last character onto the front of the string
    seq = seq[-1] + seq 
    # drop the last character (the one that was moved)
    seq = seq[0:(len(seq) - 1)] 
    # append the new sequence to the list
    bwt.append(seq) 
    # check the last character in the rotation, terminate when $
    last_character = seq[-1] 
            
# arrange list lexographically by first column
testing = bwt
bwt = sorted(bwt) ### sort this neatly

## Reformat BWT ==============================================================================

# generate count for number of rotations in bwt matrix
col_count = range(0, len(bwt))

# make empty string
bwt_last_column = ''

# iterate through each rotation in bwt matrix
for i in col_count: 
    # navigate to each rotation
    rotation = bwt[i]
    # extract the last character and assign to variable
    bwt_last_column = bwt_last_column + rotation[-1]


#print('The Burrows-Wheeler Transform is: ' + bwt_last_column)

## Reversing the BWT =====================================================================

# Generate the first column by duplicating the last, and ordering lexographically
bwt_first_column = ''.join(sorted(bwt_last_column))
    
    
# Define function for labelling sequence characters by occurance ======================================
def occurance(input):
    
    # find length of the input
    length = range(0, len(input))
    
    # set all to uppercase 
    input = input.upper()
    
    # make empty list for occurances
    order = [None] * len(input)
    
    # initialise empty counters for each base (and for special character)
    a = 0
    c = 0
    g = 0
    t = 0
    sp = 0
    
    # iterate through the input sequence
    for i in length:
        
        # check for 'A's
        if input[i] == 'A':
            a = a + 1
            order[i] = a
        
        # check for 'C's
        if input[i] == 'C':
            c = c + 1
            order[i] = c
        
        # check for 'G's
        if input[i] == 'G':
            g = g + 1
            order[i] = g
        
        # check for 'T's
        if input[i] == 'T':
            t = t + 1
            order[i] = t
        
        # check for special character '$'
        if input[i] == '$':
            sp = sp + 1
            order[i] = sp
            
    # return list containing occurances of bases in sequence          
    return order


# Concatenate first/last sequences of bwt with their base occurances
first = [str(x) + str(y) for x, y in zip(bwt_first_column, occurance(bwt_first_column))]
last = [str(x) + str(y) for x, y in zip(bwt_last_column, occurance(bwt_last_column))]

# find length of the sequence
length = range(0, len(bwt_first_column))

# Intialise the character being searched for
char_query = '$1'

# Initialise empty string for original sequence
original_sequence = ''         

# iterate through the sequence
for j in length:
    # iterate through each character/occurance pair
    for i in length:
        # check if the character in the first column matches the query
        if first[i] == char_query:
            # add the corresponding character in the last row of BWT matrix to the original sequence
            original_sequence = last[i] + original_sequence
            # update the character query with the new last column value
            char_query = last[i]
            # break from search back to main loop
            break
            

# Remove the occurance markers...
original_sequence = re.sub(r'\d+', '', original_sequence)

#... and drop the special character (now located at start position)
original_sequence = original_sequence[1:]



## Search for specific sequence within BWT ================================================

# Define search sequence
search = 'TCT'

# Set counter for number of matches
number_of_matches = 0

# Set counter for rotations
rotation_count = 0

# Set empty list for match locations
location_of_matches = []

# iterate through each bwt rotation
for rotation in testing:
    # check if sequence suffix matches to search query
    if rotation[0:len(search)] == search:
        # increase the number of matches found by 1
        number_of_matches = number_of_matches + 1
        # append the match location within the rotation to the match list
        location_of_matches.append(rotation_count)
        
    # increase rotation counter
    rotation_count = rotation_count + 1
        

print(number_of_matches)
print(location_of_matches)
print(testing) # need to use unsorted matrix, can calculate distance from $


## Make function for inverting a sequence ==================================================

def rev_complement(seq):
    
    # make new, empty list
    inverse_search = []
    
    for base in seq:
        if base == 'A':
            inverse_search.append('T')
        if base == 'G':
            inverse_search.append('C')
        if base == 'C':
            inverse_search.append('G')
        if base == 'T':
            inverse_search.append('A')
    
    # translate to str
    inverse_search = ''.join(inverse_search)
    
    # reverse str order to find reverse complement
    inverse_search = inverse_search[::-1]
    
    return inverse_search



