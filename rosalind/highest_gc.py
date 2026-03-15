def highest_gc(dnas: dict):
    highest_gc = -1
    highest_name = ''
    for name, dna in dnas.items():
        gc = gc_proportion(dna)
        if gc > highest_gc:
            highest_gc = gc
            highest_name = name
    return [highest_name, round(highest_gc, 6)]

def gc_proportion(dna):
    gc = dna.count('G') + dna.count('C')
    return (gc / len(dna)) * 100