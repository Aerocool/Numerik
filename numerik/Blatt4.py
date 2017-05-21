from math import sqrt
from pylab import *
import numpy as np
import math
from scipy import *
from numpy.linalg import cond, norm


def LU_p(A):
    A_tmp = A.copy()
    elementWithMaxValue = abs(A_tmp).max()
    n = max(A_tmp.shape)
    p = zeros((n - 1), dtype=int)

    for k in range (0, n - 1):
        i = k + argmax(abs(A_tmp[k:, k]))
        if (abs(A_tmp[i, k]) < 0.0000001 * elementWithMaxValue):
            lu = zeros(A_tmp.shape)
            return lu, p
        p[k] = i
        A_tmp[[k, i], :] = A_tmp[[i, k], :]
        A_tmp[k + 1:, k] = A_tmp[k + 1:, k] / A_tmp[k, k]
        A_tmp[k + 1:, k + 1:] = A_tmp[k + 1:, k + 1:] - outer(A_tmp[k + 1:, k], A_tmp[k, k + 1:])
    return A_tmp, p


def loeseR(U, y):
    n = len(y)
    x = zeros(n)
    x[-1] = y[-1] / U[-1, -1]
    for i in range(n - 2, -1, -1):
        x[i] = (y[i] - U[i, i + 1:].dot(x[i + 1:])) / U[i, i]
    return x


def loeseL(L, b):
    n = len(b)
    y = zeros(n, dtype=float)
    y[0] = b[0]
    for i in range(1, n):
        y[i] = b[i] - L[i, :i].dot(y[:i])
    return y



def PT(p, bb):
    b = bb.copy()
    n = len(p)
    for i in range(n - 1, -1, -1):
        b[[i, p[i]]] = b[[p[i], i]]

    return b


def loeseUT(U, y):
    n = len(y)
    x = zeros(n)
    x[0] = y[0] / U[0, 0]
    for i in range(1, n):
        x[i] = (y[i] - U[:i, i].dot(x[:i])) / U[i, i]
    return x

# Matrixmultiplikation
def mult_matrix(M, N):
    tuple_N = zip(*N)
    return [[sum(el_m * el_n for el_m, el_n in zip(row_m, col_n)) for col_n in tuple_N] for row_m in M]


def matrixMitPivot(A):
    n = len(A)

    id_mat = [[float(i ==j) for i in xrange(n)] for j in xrange(n)]

    for j in xrange(n):
        row = max(xrange(j, n), key=lambda i: abs(A[i][j]))
        if j != row:
            # Swap the rows
            id_mat[j], id_mat[row] = id_mat[row], id_mat[j]

    return id_mat


def P(p, bb):
    b = bb.copy()
    n = len(p)
    for i in range(n):
        b[[i, p[i]]] = b[[p[i], i]]

    return b


def loeseLT(L, b):
    n = len(b)
    y = zeros(n)
    y[-1] = b[-1]
    for i in range(n - 2, -1, -1):
        y[i] = b[i] - L[i + 1:, i].dot(y[i + 1:])
    return y


def lu_decomposition(A):
    n = len(A)

    L = [[0.0] * n for i in xrange(n)]
    U = [[0.0] * n for i in xrange(n)]

    P = matrixMitPivot(A)
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


def loeseQ(A, d, xx):
    x = xx.copy()
    m, n = A.shape

    for k in arange(n - 1):
        nv2 = -2.0 * A[k, k] * d[k]
        vz = A[k:, k].dot(x[k:])
        alpha = 2.0 * vz / nv2
        x[k:] = x[k:] - alpha * A[k:, k]

    return x


def loese(A, d, xx):
    x = xx.copy()
    m, n = A.shape

    for i in reversed(range(n)):
        x[i] -= A[i, i + 1:].dot(x[i + 1:])
        x[i] /= d[i]

    return x


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


def trans_matrix(M):
    n = len(M)
    return [[ M[i][j] for i in range(n)] for j in range(n)]


def normOneArgument(x):
    return sqrt(sum([x_i**2 for x_i in x]))


def Q_i(Q_min, i, j, k):
    if i < k or j < k:
        return float(i == j)
    else:
        return Q_min[i-k][j-k]


def householder(A):
    n = len(A)
    d = zeros(n, dtype=float)

    R = A
    Q = [[0.0] * n for i in xrange(n)]

    for k in range(n-1):  # We don't perform the procedure on a 1x1 matrix, so we reduce the index by 1
        I = [[float(i == j) for i in xrange(n)] for j in xrange(n)]
        x = [row[k] for row in R[k:]]
        e = [row[k] for row in I[k:]]
        alpha = -cmp(x[0],0) * normOneArgument(x)
        u = map(lambda p,q: p + alpha * q, x, e)
        norm_u = normOneArgument(u)
        v = map(lambda p: p/norm_u, u)
        Q_min = [ [float(i==j) - 2.0 * v[i] * v[j] for i in xrange(n-k)] for j in xrange(n-k) ]
        Q_t = [[ Q_i(Q_min,i,j,k) for i in xrange(n)] for j in xrange(n)]

        if k == 0:
            Q = Q_t
            R = mult_matrix(Q_t, A)
        else:
            Q = mult_matrix(Q_t, Q)
            R = mult_matrix(Q_t, R)
    return trans_matrix(Q), R, d


def hager(B):
    m, n = B.shape
    x = ones(n) / n
    iterationen = 0

    while True:
        iterationen += 1

        y = B.T.dot(x)
        xi = zeros(len(y))
        for i in range(0, len(y)-1):
            # signum funktion
            if y[i] > 0:
                xi[i] = 1
            if y[i] < 0:
                xi[i] = -1
            if y[i] == 0:
                xi[i] = 0
        z = B.dot(xi)

        if norm(z, inf) <= z.dot(x):
            return norm(y, 1)

        j = abs(z).argmax()
        x = zeros(n)
        x[j] = 1.0


def hagerinv(LU, p):
    m, n = LU.shape
    x = ones(n) / n
    iterationen = 0

    while True:
        iterationen += 1

        v = loeseUT(LU, x)
        w = loeseLT(LU, v)
        y = PT(p, w)

        xi = zeros(len(y))
        for i in range(0, len(y)-1):
            # signum funktion
            if y[i] > 0:
                xi[i] = 1
            if y[i] < 0:
                xi[i] = -1
            if y[i] == 0:
                xi[i] = 0

        v = P(p, np.array(xi))
        w = loeseL(LU, v)
        z = loeseR(LU, w)

        if norm(z, inf) <= z.dot(x):
            return norm(y, 1)

        j = abs(z).argmax()
        x = zeros(n)
        x[j] = 1.0


#Aufgabe 1

print "Aufgabe 1:"
A = array([[0.,0.,0.,1.],[2.,1.,2.,0.],[4.,4.,0.,0.],[2.,3.,1.,0.]])
LU, p = LU_p(A)

z0 = 17 + arange(A.shape[1], dtype=float)
c = A.T.dot(z0)

z = PT(p, loeseLT(LU, loeseUT(LU, c)))

print z

#Aufgabe 2

print "Aufgabe2"
for delta in 10.0 ** (-array([8, 10, 12])):
    A = array([[3, 2, 1],
               [2, 2, 2],
               [1, 2, -1]], dtype=float)
    A[1:3, 1:3] *= delta

    condinf = cond(A, inf)

    [lu, p] = LU_p(A)
    hagerCondition = hager(A) * hagerinv(lu, p)

    print delta, " ", hagerCondition, " ", condinf

#Aufgabe 5
# Matrix von Aufgabe 3
A = [[1, -5, -20], [-4, 11, -1], [8, -4, 2]]
Q, R, d = householder(A)

print "Matrix A aus Aufgabe 3:"
print A

print "Q:"
print(Q)

print "R:"
print(R)

print "Matrix aus Aufgabe 5:"
for n in [40, 50, 60]:
    A = 2.0 * eye(n) - tril(ones((n, n)))
    A[:, -1] = 1.0

    b = 2.0 - arange(n)
    b[-1] = 2.0 - n

    Q, R, d = householder(A)
    print "Q:"
    print Q
    print "R:"
    print R
    print "d:"
    print d

    print

    P, L, U = lu_decomposition(A)
    x_approx = rueckwaertseinsetzen(U, vorwaertseinsetzen(L, b))
    print x_approx
    print "x mit LU:  = ", x_approx