with open("../files/rosalind_rna.txt") as file:
    dna = file.read()
    rna = dna.replace("T", "U")
    print(rna)