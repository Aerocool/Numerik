import numpy as np
import math
from time import *
from decimal import *


a = 9.8606
c = -1.1085e25
d = 0.029


class function:
    def b(self, x):
        return a / (1.0-c*math.pow(math.e, -d*x)) - 9.0
    def bAbgeleitet(self, x):
        zaehler = -a * c * d * math.pow(math.e, d*x)
        nenner = math.pow((math.pow(math.e, d*x) - c), 2.0)
        return zaehler / nenner
    def aufgabe2(self, x):
        return x + math.log(x, math.e) - 2.0
    def aufgabe2Abgeleitet(self, x):
        return 1.0 + float(1/x)
    def aufgabe3(self, x):
        return Decimal(Decimal(math.atan(x)) - x)
    def aufgabe3Abgeleitet(self, x):
        return Decimal((1)/(Decimal(Decimal(x*x) + 1)) - 1)
    def aufgabe3AbgeleitetAbgeleitet(self, x):
        return Decimal(-(2*x/(x*x+1)**2))
    def aufgabe5X(self, x, y):
        return math.sin(x) - y
    def aufgabe5Y(self, x, y):
        return math.pow(math.e, -y) - x
    def aufgabe5XAbgeleitet(self, x, y):
        return math.cos(x)
    def aufgabe5YAbgeleitet(self, x, y):
        return -1.0*math.pow(math.e, -y)


def newtonVerfahrenMitDecimal(function, functionAbgeleitet, x_k, n):
    if n == 0:
        return x_k
    x_k = Decimal(x_k - (function(x_k)/functionAbgeleitet(x_k)))
    return newtonVerfahrenMitDecimal(function, functionAbgeleitet, x_k, n-1)


def newtonVerfahren(function, functionAbgeleitet, x_k, n):
    if n == 0:
        return x_k
    x_k = x_k - (function(x_k)/functionAbgeleitet(x_k))
    return newtonVerfahren(function, functionAbgeleitet, x_k, n-1)


def newtonVerfahren2D(functionX, functionY, functionXAbgeleitet, functionYAbgeleitet, x_k, y_k, n):
    if n == 0:
        return (x_k, y_k)
    x_kminus1 = x_k # damit bei der Iteration fuer die NST der y-Komponente nicht schon
                    # mit dem naechsten Iterationsscritt von x rechnet
    x_k = x_k - (functionX(x_k, y_k)/functionXAbgeleitet(x_k, y_k))
    y_k = y_k - (functionY(x_kminus1, y_k)/functionYAbgeleitet(x_kminus1, y_k))
    return newtonVerfahren2D(functionX, functionY, functionXAbgeleitet, functionYAbgeleitet, x_k, y_k, n-1)


def newtonVerfahrenVariante1(function, functionAbgeleitet, x_k, n, anzahlDerNullstellen):
    if n == 0:
        return x_k
    y = function(x_k)
    x_k = Decimal(y - (anzahlDerNullstellen * (function(y)/functionAbgeleitet(y))))
    return newtonVerfahrenVariante1(function, functionAbgeleitet, x_k, n-1, anzahlDerNullstellen)


def newtonVerfahrenVariante2(function, functionAbgeleitet, functionAbgeleitetAbgeleitet, x_k, n):
    if n == 0:
        return x_k
    y = function(x_k)
    x_k = Decimal(y - (function(y)*functionAbgeleitet(y))/(functionAbgeleitet(y)**2 - (function(y)*functionAbgeleitetAbgeleitet(y))))
    return newtonVerfahrenVariante2(function, functionAbgeleitet, functionAbgeleitetAbgeleitet, x_k, n-1)


def newtonVerfahrenMitFehlerabschaetzung(function, functionAbgeleitet, x_k, epsilon):
    x_kminus1 = x_k
    x_k = x_k - float(function(x_k) / functionAbgeleitet(x_k))
    if (aposteriori(0.25, x_k, x_kminus1) < epsilon):
        return x_k
    return newtonVerfahrenMitFehlerabschaetzung(function, functionAbgeleitet, x_k, epsilon)


def aposteriori(kontraktion, x_k, x_kminus1):
    return kontraktion / (1-kontraktion) * math.fabs(x_k - x_kminus1)


def sekantenVerfahren(function, x_k, x_kminus1, n):
    if n == 0:
        return x_k
    tmp = x_k
    nenner = (function(x_k) - function(x_kminus1))
    if nenner == 0: # Sobald das Sekantenverfahren nah an der NST ist, steht im Nenner eine Null, da der floating point nicht mehr genau genug ist
        return x_k
    x_k = (x_k - (x_k - x_kminus1)/nenner * function(x_k))
    x_kminus1 = tmp
    return sekantenVerfahren(function, x_k, x_kminus1, n-1)


function = function()
getcontext().prec = 70


# Aufgabe 1
timestamp = clock()
print newtonVerfahren(function.b, function.bAbgeleitet, 1961, 5)
print "Benoetigte Zeit (Newton): ", clock() - timestamp
timestamp = clock()
print sekantenVerfahren(function.b, 2000, 1961, 7)
print "Benoetigte Zeit (Sekanten): ", clock() - timestamp
# Der Zeitmessung nach bracht das Newtonverfahren minimal mehr Zeit um Maschinengenauigkeit zu erreichen,
# jedoch weniger Iterationsschritte

#Aufgabe 2c
x0 = 1
epsilon = 1e-6

nullstelle = newtonVerfahrenMitFehlerabschaetzung(function.aufgabe2, function.aufgabe2Abgeleitet, x0, epsilon)
print "Newton-Verfahren mit a-posteriori Fehlerabschaetzung: ", nullstelle
print "Probe: ", function.aufgabe2(nullstelle)

#Aufgabe 3
x0 = Decimal(1)
print "Standard-Newton_Verfahren:",\
    "\n x_1 = ", newtonVerfahren(function.aufgabe3, function.aufgabe3Abgeleitet, x0, 1),\
    "\n x_2 = ", newtonVerfahren(function.aufgabe3, function.aufgabe3Abgeleitet, x0, 2),\
    "\n x_3 = ", newtonVerfahren(function.aufgabe3, function.aufgabe3Abgeleitet, x0, 3),\
    "\n x_4 = ", newtonVerfahren(function.aufgabe3, function.aufgabe3Abgeleitet, x0, 4)

print "Variante 1:",\
    "\n x_1 = ", newtonVerfahrenVariante1(function.aufgabe3, function.aufgabe3Abgeleitet, x0, 1, Decimal(3)),\
    "\n x_2 = ", newtonVerfahrenVariante1(function.aufgabe3, function.aufgabe3Abgeleitet, x0, 2, Decimal(3)),\
    "\n x_3 = ", newtonVerfahrenVariante1(function.aufgabe3, function.aufgabe3Abgeleitet, x0, 3, Decimal(3)),\
    "\n x_4 = ", newtonVerfahrenVariante1(function.aufgabe3, function.aufgabe3Abgeleitet, x0, 4, Decimal(3))

print "Variante 2:",\
    "\n x_1 = ", newtonVerfahrenVariante2(function.aufgabe3, function.aufgabe3Abgeleitet, function.aufgabe3AbgeleitetAbgeleitet, x0, 1),\
    "\n x_2 = ", newtonVerfahrenVariante2(function.aufgabe3, function.aufgabe3Abgeleitet, function.aufgabe3AbgeleitetAbgeleitet, x0, 2),\
    "\n x_3 = ", newtonVerfahrenVariante2(function.aufgabe3, function.aufgabe3Abgeleitet, function.aufgabe3AbgeleitetAbgeleitet, x0, 3),\
    "\n x_4 = ", newtonVerfahrenVariante2(function.aufgabe3, function.aufgabe3Abgeleitet, function.aufgabe3AbgeleitetAbgeleitet, x0, 4)

#Aufgabe 5

x0 = 0
y0 = 1

print "Newton-Verfahren fuer die Aufgabe 5: (x,y) = ",\
    newtonVerfahren2D(function.aufgabe5X, function.aufgabe5XAbgeleitet,\
                      function.aufgabe5Y, function.aufgabe5YAbgeleitet, x0, y0, 10)