import manager
import hamming, highest_gc

# dnas = manager.getstr('rosalind_hamm.txt').split()
# print(hamming.hammingdistance(dnas[0], dnas[1]))

# Highest GC content
text = manager.getstr('rosalind_gc.txt')
dnas = manager.fasta_to_dna(text)
print(highest_gc.highest_gc(dnas))