def Dec2Dna(y):
    # Initialization
    n = len(y)
    dna = []

    # Converts and write into the file
    for i in range(n):
        if y[i] == 0:
            dna.append('A')
        elif y[i] == 1:
            dna.append('C')
        elif y[i] == 2:
            dna.append('G')
        elif y[i] == 3:
            dna.append('T')
        else:
            raise ValueError("Unexpected symbol, outside the quaternary alphabet (0,1,2,3)")

    return dna
