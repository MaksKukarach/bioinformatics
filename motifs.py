"""
Count(Motifs)

Return per-column counts of nucleotides for a list of equal-length DNA strings.

Args:
    Motifs (list of str): list of strings of equal length using characters 'A','C','G','T'.

Returns:
    dict: mapping each nucleotide 'A','C','G','T' to a list of integer counts for each column.

Example:
    >>> Count(["AT","AC","GT"])
    {'A': [2, 0], 'C': [0, 1], 'G': [1, 0], 'T': [0, 2]}

"""
def Count(Motifs):
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
def Profile(Motifs):
    n = len(Motifs)
    k = len(Motifs[0])
    counts = Count(Motifs)
    for lst in counts.values():
        for i in range(k):
            lst[i] /= n
            lst[i] = round(lst[i], 2)
    return counts