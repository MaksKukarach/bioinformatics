import random

"""
Count(Motifs: list[str])

Return per-column counts of nucleotides for a list of equal-length DNA strings.

Args:
    Motifs (list of str): list of strings of equal length using characters 'A','C','G','T'.

Returns:
    dict: mapping each nucleotide 'A','C','G','T' to a list of integer counts for each column.

Example:
    >>> Count(["AT","AC","GT"])
    {'A': [2, 0], 'C': [0, 1], 'G': [1, 0], 'T': [0, 2]}

"""
def Count(Motifs: list[str]):
    # Find rows (n) and colums (k)
    n = len(Motifs)
    k = len(Motifs[0])
    # Set up the count dictionary
    count = {}
    for nucleotide in 'ACGT':
        count[nucleotide] = [0] * k
    # Count
    for i in range(n):
        for j in range(k):
            nucleotide = Motifs[i][j]
            count[nucleotide][j] += 1
    return count

"""
Profile(Motifs: list[str])

Compute the profile (column-wise nucleotide frequencies) for a set of DNA motifs.

Args:
    Motifs (list of str): equal-length DNA strings using characters 'A', 'C', 'G', 'T'.

Returns:
    dict: mapping each nucleotide ('A','C','G','T') to a list of floats giving the frequency
    of that nucleotide in each column.
Example:
    >>> Profile(["AT", "AC", "GT"])
    {'A': [0.67, 0.0], 'C': [0.0, 0.33], 'G': [0.33, 0.0], 'T': [0.0, 0.67]}
"""
def Profile(Motifs: list[str], pseudocounts=False):
    if pseudocounts:
        return ProfileWithPseudocounts(Motifs)
    n = len(Motifs)
    k = len(Motifs[0])
    counts = Count(Motifs)
    for lst in counts.values():
        for i in range(k):
            lst[i] /= n
    return counts

"""
Consensus(Motifs: list[str])

Return the consensus string for a list of equal-length DNA motifs.

Args:
    Motifs (list[str]): DNA strings of equal length using 'A', 'C', 'G', 'T'.
Returns:
    str: consensus sequence (most frequent nucleotide per column).
Example:
    >>> Consensus(["AT", "AC", "GT"])
    'AT'
"""
def Consensus(Motifs: list[str]):
    count = Count(Motifs)
    k = len(Motifs[0])
    consensus = ""
    for i in range(k):
        m = 0
        mostFrequentSymbol = ""
        for symbol in 'ACGT':
            if count[symbol][i] > m:
                m = count[symbol][i]
                mostFrequentSymbol = symbol
        consensus += mostFrequentSymbol
    return consensus

def Score(Motifs: list[str]) -> int:
    consensus = Consensus(Motifs)
    count = Count(Motifs)
    k = len(Motifs[0])
    score = 0
    for i in range(k):
        for symbol in 'ACGT':
            if symbol != consensus[i]:
                score += count[symbol][i]
    return score

def Probability(Sequence, Profile):
    probability = 1
    for i in range(len(Sequence)):
        symbol = Sequence[i]
        probability *= Profile[symbol][i]
    return probability

def ProfileMostProbableKmer(text, k, profile):
    most_probable = ""
    highest_probability = -1
    for i in range(0, len(text) - k + 1):
        sequence = text[i:i+k]
        probability = Probability(sequence, profile)
        if probability > highest_probability:
            highest_probability = probability
            most_probable = sequence
    return most_probable

def GreedyMotifSearch(DNAs, k, t, pseudocounts=False):
    BestMotifs = []
    # Initial best guess
    for i in range(t):
        BestMotifs.append(DNAs[i][0:k])
    BestScore = Score(BestMotifs)
    for i in range(len(DNAs[0]) - k + 1):
        # Take a k-mer in DNAs[0]
        Motifs = [DNAs[0][i:i+k]]
        # For each next DNA string, find a k-mer that is 
        # Profile-most probable for our current profile matrix.
        for j in range(1, t):
            ProfileMatrix = Profile(Motifs, pseudocounts)
            kmer = ProfileMostProbableKmer(DNAs[j], k, ProfileMatrix)
            Motifs.append(kmer)
        # After we picked our list of best k-mers, calculate its score
        # and compare to the current BestMotifs
        NewScore = Score(Motifs)
        if NewScore < BestScore: # Important: the lower the score, the better!
            BestScore = NewScore
            BestMotifs = Motifs
    return BestMotifs

def CountWithPseudocounts(Motifs):
    counts = Count(Motifs)
    for lst in counts.values():
        for i in range(len(lst)):
            lst[i] += 1
    return counts

def ProfileWithPseudocounts(Motifs):
    k = len(Motifs[0])
    counts = CountWithPseudocounts(Motifs)
    n = len(Motifs) + 4 # Note that since we added 1 to each letter, we must add 4 to total!
    for lst in counts.values():
        for i in range(k):
            lst[i] /= n
    return counts

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    return GreedyMotifSearch(Dna, k, t, pseudocounts=True)

# Get a list with one profile-most probable kmer from each Dna string
def Motifs(Profile, DNAs):
    Motifs = []
    k = len(list(Profile.values())[0])
    for Dna in DNAs:
        Motifs.append(ProfileMostProbableKmer(Dna, k, Profile))
    return Motifs

def RandomMotifs(DNAs, k):
    Motifs = []
    for Dna in DNAs:
        start = random.randint(0, len(Dna) - k)
        motif = Dna[start:start+k]
        Motifs.append(motif)
    return Motifs

def RandomizedMotifSearch(DNAs, k):
    best_motifs = RandomMotifs(DNAs, k)
    while True:
        profile = ProfileWithPseudocounts(best_motifs)
        motifs = Motifs(profile, DNAs)
        if Score(motifs) < Score(best_motifs): # Note: the score should be lower for it to be better
            best_motifs = motifs
        else:
            return best_motifs

# Normalize numbers in a dict so that sum up to 1
def Normalize(Probabilities: dict[str, float]) -> dict[str, float]:
    Sum = sum(Probabilities.values())
    for key in Probabilities:
        # Divide each number by the sum of all numbers
        Probabilities[key] /= Sum
    return Probabilities

# My original function
# Assumes the sum of Probabilities is 1
def DeprecatedWeightedDie(Probabilities):
    target = random.uniform(0, 1)
    Probabilities = sorted(list(Probabilities.items()), key=lambda pair: pair[1])
    lower_bound = 0
    for item, probability in Probabilities:
        upper_bound = lower_bound + probability
        if target >= lower_bound and target <= upper_bound:
            return item
        lower_bound = upper_bound
    return None

# Updated function
# Assumes the sum of Probabilities is 1
def WeightedDie(Probabilities):
    p = random.uniform(0, 1)
    upperbound = 0
    for key, prob in Probabilities.items():
        upperbound += prob
        if p < upperbound:
            return key
    return None