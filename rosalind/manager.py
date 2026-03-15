def getstr(filename):
    with open(f'../files/{filename}') as file:
        str = file.read().strip()
        return str

def fasta_to_dna(text):
    items = text.split('>')
    dnas = dict()
    for item in items:
        item = item.split('\n')
        name = item[0]
        dna = ''.join(item[1:])
        dnas[name] = dna
    return dnas

