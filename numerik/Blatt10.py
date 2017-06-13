from sympy import *
import numpy as np


class function:
    def aufgabe1(self, x):
        return 1.0/(1.0+x*x)
    def aufgabe5(self, x):
        return np.sin(np.pi*(x*x))


def summierteTrapezregel(m, b, a, f):
    l = float(b - a) / float(m)
    a_schlange = [None] * (m + 1)
    result = 0.0
    for i in range(0, m + 1):
        a_schlange[i] = a + i * l
    for i in range(0, len(a_schlange)):
        if i == 0 or i == len(a_schlange) - 1:
            result += f(a_schlange[i])
        else:
            result += 2 * f(a_schlange[i])

    return result * (l / 2.0)


def summierteSimpsonRegel(m, b, a, f):
    l = float(b-a)/float(m)
    a_schlange = [None] * (m+1)
    result = 0.0
    for i in range(0, m+1):
        a_schlange[i] = a + i*l
    for i in range(0, len(a_schlange)):
        if i == 0 or i == len(a_schlange)-1:
            result += f(a_schlange[i])
        elif (i % 2) == 0:
            result += 2 * f(a_schlange[i])
        else:
            result += 4 * f(a_schlange[i])
    return result * (l/3.0)


def gaussVerfahren(betas, xs, f):
    result = 0.0
    for i in range(0, len(betas)):
        result += betas[i] * f(xs[i])
    return result


def romberg(f, a, b, n):
    r = np.array( [[0] * (n+1)] * (n+1), float )
    h = b - a
    r[0,0] = 0.5 * h * ( f( a ) + f( b ) )
    powerOf2 = 1
    for i in xrange( 1, n + 1 ):
        h = 0.5 * h
        sum = 0.0
        powerOf2 = 2 * powerOf2
        for k in xrange( 1, powerOf2, 2 ):
            sum = sum + f( a + k * h )
        r[i,0] = 0.5 * r[i-1,0] + sum * h
        powerOf4 = 1
        for j in xrange( 1, i + 1 ):
            powerOf4 = 4 * powerOf4
            r[i,j] = r[i,j-1] + ( r[i,j-1] - r[i-1,j-1] ) / ( powerOf4 - 1 )
    return r


def bulirsch(F, x, y, xStop, tol):
    def midpoint(F, x, y, xStop, nSteps):
        h = (xStop - x)/ nSteps
        y0 = y
        y1 = y0 + h*F(x, y0)
        for i in range(nSteps-1):
            x = x + h
            y2 = y0 + 2.0*h*F(x, y1)
            y0 = y1
            y1 = y2
        return 0.5*(y1 + y0 + h*F(x, y2))


#Aufgabe 1
# Teil a
x = symbols("x")
f = 1/(1+x**2) #Funktion
print f.integrate(x)
F = Lambda(x, f.integrate(x))
print "Exaktes Ergebnis:", F(1) - F(0)
ergebnis = np.pi/4.0
print "Exaktes Ergebnis als float:", ergebnis

f = function()
# Teil c
ergebnisSimpson = summierteSimpsonRegel(4, 1, 0, f.aufgabe1)
print "Ergebnis mit der summierten Simpson-Regel:", ergebnisSimpson
print "Absoluter Fehler mit der summierten Simpson-Regel:", np.abs(ergebnisSimpson - ergebnis)
# Teil b
ergebnisTrapez = summierteTrapezregel(8, 1, 0, f.aufgabe1)
print "Ergebnis mit der summierten Trapez-Regel:", ergebnisTrapez
print "Absoluter Fehler mit der summierten Trapez-Regel:", np.abs(ergebnisTrapez - ergebnis)
# Teil d
xs = [-np.math.sqrt(3.0/5.0), 0, np.math.sqrt(3.0/5.0)]
betas = [5.0/9.0, 8.0/9.0, 5.0/9.0]
ergebnisGauss = gaussVerfahren(betas, xs, f.aufgabe1)
print "Ergebnis mit Gauss:", ergebnisGauss
print "Absoluter Fehler mit Gauss:", np.abs(ergebnisGauss - ergebnis)

#Aufgabe 5
m = [1, 2, 3, 4]
for i in range(0, len(m)):
    print "Trapezregel mit m = ", m[i], "ist:", summierteTrapezregel(m[i], 1, -1, f.aufgabe5)
ergebnisRomberg = romberg(f.aufgabe5, -1.0, 1.0, 4) # Hier sind alle Zwischenwerte
print "Mit Romberg:", ergebnisRomberg[len(ergebnisRomberg)-1][len(ergebnisRomberg)-1]
#print "Mit Bulrisch:", Funktion ist oben implementiert

#Aufgabe 1b Skript S.223

