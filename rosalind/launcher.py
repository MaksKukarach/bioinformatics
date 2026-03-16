import manager
import hamming, highest_gc, rna_to_protein, dna

# dnas = manager.getstr('rosalind_hamm.txt').split()
# print(hamming.hammingdistance(dnas[0], dnas[1]))

# Highest GC content
# text = manager.getstr('rosalind_gc.txt')
# dnas = manager.fasta_to_dna(text)
# print(highest_gc.highest_gc(dnas))

# RNA -> Protein
# rna = manager.getstr('rosalind_prot.txt')
# protein = rna_to_protein.rna_to_protein(rna)
# print(protein)

# Motif (substring) positions
# my_dna, my_motif = manager.getstr('rosalind_subs.txt').split('\n')
# print(my_dna, my_motif)
# print(*dna.motif_positions(my_dna, my_motif))