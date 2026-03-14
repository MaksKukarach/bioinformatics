with open("../files/rosalind_revc.txt") as file:
    dna = file.read().strip()
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    comp = ""
    for letter in dna:
        comp += pairs[letter]
    comp = comp[::-1]
    print(comp)
