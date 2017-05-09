import numpy as np
import math


def getResiduum(b, A, x_approx):
    r = b - A.dot(x_approx)
    for i in range(0, len(r)):
        r[i] = abs(r[i])
    return r


def getS(b, A, x_approx):
    b = abs(b)
    A = abs(A)
    x_approx = abs(x_approx)
    s = A.dot(x_approx) + b
    return s


def calculateEpsilon(r, s):
    epsilon = float(1)
    for i in range(0, len(r)):
        if r[i] > s[i]:
            tmp = float(r[i] / s[i])
            if tmp > epsilon:
                epsilon = tmp
    return epsilon


def mult_matrix(M, N):
    tuple_N = zip(*N)
    return [[sum(el_m * el_n for el_m, el_n in zip(row_m, col_n)) for col_n in tuple_N] for row_m in M]


def pivot_matrix(M):
    m = len(M)

    id_mat = [[float(i ==j) for i in xrange(m)] for j in xrange(m)]

    for j in xrange(m):
        row = max(xrange(j, m), key=lambda i: abs(M[i][j]))
        if j != row:
            # Swap the rows
            id_mat[j], id_mat[row] = id_mat[row], id_mat[j]

    return id_mat


def lu_decomposition(A):
    n = len(A)

    L = [[0.0] * n for i in xrange(n)]
    U = [[0.0] * n for i in xrange(n)]

    P = pivot_matrix(A)
    PA = mult_matrix(P, A)

    for j in xrange(n):
        L[j][j] = 1.0

        for i in xrange(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in xrange(i))
            U[i][j] = PA[i][j] - s1

        for i in xrange(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in xrange(j))
            L[i][j] = (PA[i][j] - s2) / U[j][j]

    return (P, L, U)


def vorwaertseinsetzen(L, b):
    y = [0.0] * len(b)
    for i in range(0, len(L[0])):
        y[i] = b[i]
        for k in range(0, i):
            tmp = y[k] * L[i][k]
            y[i] = y[i] - tmp
    return y


def rueckwaertseinsetzen(U, y):
    x = [0.0] * len(y)
    x[len(U[0])-1] = float(y[len(y)-1] / U[len(U[0])-1][len(U[0])-1])
    for i in range(len(U[0])-2, -1, -1):
        tmp = 0.0
        for k in range (len(U[0])-1, -1, -1):
            if k <= i:
                break
            tmp = tmp + U[i][k] * x[k]
        x[i] = float((y[i] - tmp) / U[i][i])
    return x


# Aufgabe4: Gleichungssystem mit Matrix A aus Aufgabe 2 loesen:
beta = 10
n = [10, 15, 20]

for i in n:
    A = [[0]*i for _ in range(i)]
    for k in range(0, i):
        for j in range(0, i):
            if k == j and k <> (i-1):
                A[k][j] = 1
                A[k+1][j] = -beta
            if k == 0 and j == (i-1):
                A[k][j] = beta
    b = np.zeros(i)
    for k in range(1, i):
        b[k] = 1 - beta
    b[0] = 1 + beta
    b[i-1] = -beta

    P, L, U = lu_decomposition(A)
    print "Loesung: ", rueckwaertseinsetzen(U, vorwaertseinsetzen(L, b))

# Aufgabe5
delta = [1*math.pow(10., -8.), 1*math.pow(10., -10.), 1*math.pow(10., -12.)]
for i in delta:
    A = [[3, 2, 1], [2, 2*i, 2*i], [1, 2*i, -i]]
    b = [3 + 3*i, 6*i, 2*i]
    P, L, U = lu_decomposition(A)
    x_approx = rueckwaertseinsetzen(U, vorwaertseinsetzen(L, b))
    print x_approx
    x_approx = np.array(x_approx)
    A = np.array(A)
    b = np.array(b)

    print "Fuer delta = ", i
    print "|r| = ", getResiduum(b, A, x_approx)
    print "s = ", getS(b, A, x_approx)
    print "epsilon: ", calculateEpsilon(getResiduum(b, A, x_approx), getS(b, A, x_approx))
    print "cond(A): ", np.linalg.cond(A)

    # Die Approximation von x kann bei jedem Delta mit dem Verfahren von Prager und Oettli akzeptiert werden
