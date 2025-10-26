def MinimumSkew(Genome):
    positions = [] # output variable
    skew = SkewArray(Genome)
    minimum = min(skew)
    for i in range(len(skew)):
        if skew[i] == minimum:
            positions.append(i)
    return positions


# Input:  A String Genome
# Output: The skew array of Genome as a list.
def SkewArray(Genome):
    skew = [0]
    for i in range(len(Genome)):
        if Genome[i] == "C":
            skew.append(skew[i] - 1)
        elif Genome[i] == "G":
            skew.append(skew[i] + 1)
        else:
            skew.append(skew[i])
    return skew


# Create a window and slide it, updating counts according to first and next symbols
def FasterSymbolMap(Genome, symbol):
    mapping = {}

    # 1.Extend the genome to account for circularity
    ExtendedGenome = Genome + Genome[:len(Genome)//2]
    # 2.Create first window count
    len_window = len(Genome)//2
    mapping[0] = PatternCount(ExtendedGenome[:len_window], symbol)
    # 3.Slide and update last and next symbols
    for i in range(1, len(Genome)):
        count = mapping[i-1]
        if ExtendedGenome[i-1] == symbol:
            count -= 1
        if ExtendedGenome[i+len_window-1] == symbol:
            count += 1
        mapping[i] = count

    return mapping

def SymbolMap(Genome, symbol):
    mapping = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        mapping[i] = PatternCount(ExtendedGenome[i:i+(n//2)], symbol)
    return mapping

def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 

# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words
# Copy your FrequencyMap() function here.
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] += 1
    return freq

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