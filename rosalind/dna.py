def dna_to_rna(dna):
    rna = dna.replace("T", "U")
    return rna

def reverse_comp(dna):
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    comp = ""
    for letter in dna:
        comp += pairs[letter]
    return comp[::-1]

def motif_positions(dna: str, motif: str, one_based_indexing=True):
    """
    Return motif positions, using one-based indexing by default.
    """
    positions = []
    for i in range(0, len(dna) - len(motif) + 1):
        if dna[i:i+len(motif)] == motif:
            positions.append(i)
    if one_based_indexing:
        positions = [i+1 for i in positions]
    return positions