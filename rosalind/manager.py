import hamming

def getstr(filename):
    with open(f'../files/{filename}') as file:
        str = file.read().strip()
        return str
    
dnas = getstr('rosalind_hamm.txt').split()
print(hamming.hammingdistance(dnas[0], dnas[1]))