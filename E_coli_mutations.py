import replication

with open("files/E_coli.txt", "r") as file:
    genome = file.read()

C_counts = replication.SymbolArray(genome, "C")
print(C_counts)