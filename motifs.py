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
    of that nucleotide in each column, rounded to 2 digits.
Example:
    >>> Profile(["AT", "AC", "GT"])
    {'A': [0.67, 0.0], 'C': [0.0, 0.33], 'G': [0.33, 0.0], 'T': [0.0, 0.67]}
"""
def Profile(Motifs: list[str]):
    n = len(Motifs)
    k = len(Motifs[0])
    counts = Count(Motifs)
    for lst in counts.values():
        for i in range(k):
            lst[i] /= n
            lst[i] = round(lst[i], 2)
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