from sympy import  *


def summierteTrapezregel():
    return 0


def summierteSimpsonRegel():
    return 0


def gaussVerfahren():
    return 0


#Aufgabe 1a
x = symbols("x")
f = 1/(1+x**2) #Funktion
print f.integrate(x)
F = Lambda(x, f.integrate(x))
print "Exaktes Ergebnis:", F(1) - F(0)


#Aufgabe 1b Skript S.223

