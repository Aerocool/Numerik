import numpy as np

def zerlegung(A):
    for spalte in range(0, len(A[0])):
        if spalte == len(A[0]) - 1:
            continue
        if A[spalte][spalte] == 0:
            for r in range(spalte, len(A)):
                if A[r][spalte] <> 0:
                    A = vertauscheZeilen(A, r, spalte)
        l = spalte + 1
        for l in range(0, len(A[0])):
            for i in range(0, len(A[l])):
                A[l][i] = float(A[l][i] / A[spalte][spalte])
    return A

def permutation(p, x):
    return 0

def vertauscheZeilen(A, z1, z2):
    P = np.zeros((len(A), len(A))) # P = Permutationsmatrix
    for i in range (0, len(P)): # Zeilen
        for j in range (0, len(P)): # Spalten
            if (i == z1) and (j == z2):
                P[i][j] = 1
            elif (i == z2) and (j == z1):
                P[i][j] = 1
            elif (i == j) and (i <> z1) and (i <> z2):
                P[i][j] = 1
    return np.mat(P) * np.mat(A)

def vorwaerts(LU, x):
    y = np.zeros(len(x))
    for i in range (0, len(y)):
        if i > 0:
            for j in range(0, i):
                y[i] += LU[i][j] * y[j]
        else:
            y[i] = x[i]
    return y

def ruecwaerts(LU, x):
    return 0

#l = np.array(((2, 3), (3, 5)))
#a = np.array(((1, 2), (5, -1)))

#print np.dot(l, a)
#print np.mat(l) * np.mat(a)

A = np.array(((2, 1, 1), (4, 4, 0), (2, 7, 1)))
#result = zerlegung(A)
result = vertauscheZeilen(A, 0, 1)
print result
