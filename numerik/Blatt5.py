import numpy as np
from pylab import *


def Gauss_Seidel(A, b, x_kplus1, steps, richtigeLoesung):
    i = 0
    j = 0
    if steps == 0:
        return x_kplus1

    for i in range(0, len(b)):
        x_k = np.zeros(len(x_kplus1))
        for l in range (0, len(x_kplus1)):
            x_k[l] = x_kplus1[l]
        ersteSumme = 0
        zweiteSumme = 0
        for j in range(0, i):
            ersteSumme += A[i][j] * x_kplus1[j]
        for j in range(i+1, len(A)):
            zweiteSumme += A[i][j] * x_k[j]
        x_kplus1[i] = float(1.0 / A[i][i]) * (b[i] - ersteSumme - zweiteSumme)

    fehler = np.zeros(len(x_kplus1))
    for i in range(0, len(x_kplus1)):
        fehler[i] = norm(x_kplus1[i] - richtigeLoesung[i])
    print "|| x_k - x || = ", fehler
    return Gauss_Seidel(A, b, x_kplus1, steps - 1, richtigeLoesung)


#Aufgabe 2
A = [[2, 0, 1], [0, 1, 0], [1, 0, 2]]
b = [4, 0, 5]
x = [0, 0, 0]
richtigeLoesung = [1, 0, 2]

x = Gauss_Seidel(A, b, x, 10, richtigeLoesung)
print "Ergebnis aus Aufgabe 1 als Testbeispiel: ", x

A = np.zeros((30,) * 2)
b = np.zeros(30)
richtigeLoesung = np.ones(30)


for i in range (0, len(A)):
    A[i][i] = 1
    if i <> (len(A)-1):
        A[i][i+1] = float(-1.16)
        A[i+1][i] = float(0.16)

b[0] = float(-0.16)
b[29] = float(1.16)
x_null = np.zeros(30)

x_null = Gauss_Seidel(A, b, x_null, 120, richtigeLoesung)

print "Ergebnis: ", x_null

