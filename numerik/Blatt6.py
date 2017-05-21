from pylab import *
from numpy import *
import numpy as np


def squareSum(A):
    retturnValue = 0
    for i in range(A.shape[0]):
        for j in range(i, A.shape[1], 1):
            if i != j:
                retturnValue += A[i][j] * A[i][j]
    return retturnValue


def sign(a):
    if (a >= 0):
        return 1
    else:
        return -1



def vectorIt(A, n, x0):
    x = x0
    for i in range(0, n, 1):
        y = A.dot(x)
        xn = y / norm(y)
        p = x.dot(y)
        x = xn
    return p


def rotationAlaJacobi(A, precision):
    qSum = precision + 1
    Qi = eye(A.shape[0], dtype=float)
    while qSum > precision:
        upmat = abs(triu(A, 1))
        m = argmax(upmat)
        mIndex = unravel_index(m, A.shape)
        i = mIndex[0]
        j = mIndex[1]
        aij = A[i][j]
        aii = A[i][i]
        alpha = (aij - aii) / (2. * aij)
        c = math.sqrt(1 / 2. + ((1 / 2.) * ((alpha * alpha) / (1. + alpha * alpha))))
        s = sign(alpha) / (2. * c * math.sqrt(1 + alpha * alpha))
        Qi = Qi.dot(rotationAlaGivens(A.shape[0], i, j, c, s))
        A = ((Qi.transpose()).dot(A)).dot(Qi)
        qSum = squareSum(triu(A, 1))
    return A


def rotationAlaGivens(size, i0, j0, c, s):
    Q = eye(size, dtype=float)
    Q[i0][j0] = s
    Q[j0][i0] = -s
    Q[i0][i0] = c
    Q[j0][j0] = c
    return Q



from numpy.linalg import qr


def QRIt(A, epsilon):
    Ak = A
    qSum = epsilon + 1
    while qSum > epsilon:
        Q, R = qr(Ak)
        Ak = R.dot(Q)
        qSum = squareSum(triu(Ak, 1))
    return Ak


def getMatrixB(size):
    B = eye(size, dtype=float)
    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            if (i == j):
                B[i][j] = 4
                i += B.shape[0] - 1
                break
            if (i == j + 1):
                B[i][j] = -1
                B[j][i] = -1
                j += B.shape[0] - 1
    return B


def getAPartOfMatrix(B):
    I = eye(B.shape[0], dtype=float)
    returnMatrix = eye(B.shape[0] * B.shape[0], dtype=float)
    for i in range(returnMatrix.shape[0]):
        for j in range(returnMatrix.shape[1]):
            if ((i == j) and ((j % (B.shape[0])) == 0)):
                for k in range(B.shape[0]):
                    for j in range(B.shape[1]):
                        returnMatrix[i + k][j + j] = B[k][j]
                i += B.shape[0] - 1
                break
            if ((i == j + B.shape[0]) and ((i % (B.shape[0])) == 0)):
                for k in range(B.shape[0]):
                    for j in range(B.shape[1]):
                        returnMatrix[i + k][j + j] = (-1) * I[k][j] + 0
                        returnMatrix[j + j][i + k] = (-1) * I[k][j] + 0
                j += B.shape[0] - 1
    return returnMatrix


#Aufgabe 1
A = eye(3, dtype=float)
A[0,0]=4
A[0,1]=2
A[0,2]=1
A[1,0]=2
A[1,1]=4
A[1,2]=2
A[2,0]=1
A[2,1]=2
A[2,2]=4

n = 5
x0 = zeros(3)
x0[0] = 1
result1 = vectorIt(A, n, x0)
print "Aufgabe 1:"
print result1

#Aufgabe 2
epsilon = 10 ** (-3)
B = getMatrixB(10)
A = getAPartOfMatrix(B)
rotation = rotationAlaJacobi(A, epsilon)
result2 = [rotation[i][i] for i in range(len(rotation))]
result2.sort()
print "Aufgabe 2:"
print result2

#Aufgabe 4
epsilon = 10 ** (-3)
B = getMatrixB(10)
A = getAPartOfMatrix(B)
result3 = QRIt(A, epsilon)
print "Aufgabe 3:"
print result3
