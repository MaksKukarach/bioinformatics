from manager import getstr, fasta_to_dna
import hamming, highest_gc, rna_to_protein, dna, heredity

# dnas = getstr('rosalind_hamm.txt').split()
# print(hamming.hammingdistance(dnas[0], dnas[1]))

# Highest GC content
# text = getstr('rosalind_gc.txt')
# dnas = fasta_to_dna(text)
# print(highest_gc.highest_gc(dnas))

# RNA -> Protein
# rna = getstr('rosalind_prot.txt')
# protein = rna_to_protein.rna_to_protein(rna)
# print(protein)

# Motif (substring) positions
# my_dna, my_motif = getstr('rosalind_subs.txt').split('\n')
# print(my_dna, my_motif)
# print(*dna.motif_positions(my_dna, my_motif))

# # DNA palindromes
# text = getstr('rosalind_revp.txt')
# # text = '>Rosalind_24\nTCAATGCATGCGGGTCTATATGCAT'
# my_dna = list(fasta_to_dna(text).values())[0]
# pals = dna.reverse_palindromes(my_dna)
# for pair in pals:
#     print(pair[0], pair[1])

# Probability of dominant allele given a sample
# print(heredity.prob_dominant(26, 25, 17))