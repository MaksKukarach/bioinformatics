import manager
import hamming

# dnas = manager.getstr('rosalind_hamm.txt').split()
# print(hamming.hammingdistance(dnas[0], dnas[1]))

text = manager.getstr('rosalind_gc.txt')
dnas = manager.fasta_to_dna(text)
print(dnas)