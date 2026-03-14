with open("../files/rosalind_dna.txt") as f:
    dna = f.read()
    print(dna)
    print(dna.count("A"), dna.count("C"), dna.count("G"), dna.count("T"))