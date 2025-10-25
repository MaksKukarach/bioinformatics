import replication
import json

# --- Read genome ---
with open("files/E_coli.txt", "r") as file:
    genome = file.read().strip()

# --- Compute symbol map ---
count_mapping = replication.FasterSymbolMap(genome, "C")

# --- Save to JSON ---
with open("files/count_mapping.json", "w") as file:
    json.dump(count_mapping, file, indent=4)