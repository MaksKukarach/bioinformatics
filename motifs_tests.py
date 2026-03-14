from motifs import *

##### motifs.py #####
my_motifs = ["TATA", "TATC", "GATA", "TTTA", "TAAA"]
print(Count(my_motifs))
print(
    Profile(my_motifs)
)  # Yes, they do sum up to 1, if you look at sum of i-th elements in lists,
# not the sum of values in each list!

odd_frequencies = ["AT", "AC", "GT"]
print("Odd frequencies: ")
print(Profile(odd_frequencies))

profile = {
    "A": [0.2, 0.2, 0.3, 0.2, 0.3],
    "C": [0.4, 0.3, 0.1, 0.5, 0.1],
    "G": [0.3, 0.3, 0.5, 0.2, 0.4],
    "T": [0.1, 0.2, 0.1, 0.1, 0.2],
}
print(
    ProfileMostProbableKmer(
        "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", 5, profile
    )
)

DNAs = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
k = 3
t = 5
result = GreedyMotifSearch(DNAs, k, t)
expected = ["CAG", "CAG", "CAA", "CAA", "CAA"]
assert result == expected, f"GreedyMotifSearch got {result}, expected {expected}"
print("GreedyMotifSearch passed!")
