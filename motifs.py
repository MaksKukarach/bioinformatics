# Input:  A list of n strings, each representing a k-mer (Motif)
# Output: Count of nucleotides by position in k-mers
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

def Profile(Motifs):
    n = len(Motifs)
    k = len(Motifs[0])
    counts = Count(Motifs)
    for lst in counts.values():
        for i in range(k):
            lst[i] /= n
    return counts