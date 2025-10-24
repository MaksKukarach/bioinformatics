# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(Pattern):   
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

# Copy your Reverse() function here.
def Reverse(Pattern):
    reverse = ""
    for letter in Pattern:
        reverse = letter + reverse
    return reverse

# Copy your Complement() function here.
def Complement(Pattern):
    complements = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    complement = ''
    for letter in Pattern:
        if letter in complements.keys():
            complement += complements[letter]
    return complement