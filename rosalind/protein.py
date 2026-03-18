codon_table = {
    "UUU": "F", "UUC": "F", 
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", 
    "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", 
    "CAU": "H", "CAC": "H",
    "AAU": "N", "AAC": "N", 
    "GAU": "D", "GAC": "D",
    "UAA": "Stop", "UAG": "Stop", 
    "CAA": "Q", "CAG": "Q",
    "AAA": "K", "AAG": "K", 
    "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", 
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", 
    "AGU": "S", "AGC": "S",
    "AGA": "R", "AGG": "R", 
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", 
    "UGA": "Stop", 
    "UGG": "W"
}

def rna_to_protein(rna):
    result = ""
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        result += codon_table[codon]
    return result

aa_weights = {
    "A": 71.03711,
    "C": 103.00919,
    "D": 115.02694,
    "E": 129.04259,
    "F": 147.06841,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "K": 128.09496,
    "L": 113.08406,
    "M": 131.04049,
    "N": 114.04293,
    "P": 97.05276,
    "Q": 128.05858,
    "R": 156.10111,
    "S": 87.03203,
    "T": 101.04768,
    "V": 99.06841,
    "W": 186.07931,
    "Y": 163.06333
}

def weight(protein):
    weight = 0
    for aa in protein:
        weight += aa_weights[aa]
    return weight