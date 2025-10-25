import json
import matplotlib.pyplot as plt

with open("files/count_mapping.json", "r") as file:
    count_mapping = json.load(file)

# --- Prepare data for plotting ---
positions = list(map(int, count_mapping.keys()))
counts = list(map(int, count_mapping.values()))

# --- Plot ---
plt.figure(figsize=(10, 5))
plt.plot(positions, counts, color="blue", linewidth=1)
plt.title("C counts across E. coli genome")
plt.xlabel("Position")
plt.ylabel("Count of 'C' in half-genome window")
plt.grid(True)
plt.tight_layout()
plt.show()