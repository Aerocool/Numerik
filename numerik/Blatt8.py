import numpy as np
import math
import matplotlib.pyplot as plt
import pylab

class function:
    def aufgabe2b(self, x):
        return float(1.0 / 1.0 + x*x)
    def aufgabe2c(self, x):
        return float(1/(1+x*x))
    def aufgabe2dX(self, t):
        return (1.0/(1+t)) * math.cos(3*math.pi*t)
    def aufgabe2dY(self, t):
        return (1.0/(1+t)) * math.sin(3*math.pi*t)
    def mAufgabe2b(self, i, m):
        return float(-5.0 + float(10/(m-1))*i)
    def mAufgabe2c(self, i, m):
        return float(-5.0 * math.cos(math.pi*(2.0*i+1)/(2*m)))
    def mAufgabe2d(self, i, m):
        return float(i/(m-1.0))

def dividierteDifferenzen (y, x):
    d = np.empty((len(y), len(y)))
    for i in range(0, len(d)):
        for j in range(0, len(d)):
            if i == j:
                d[i][i] = y[i]
            else:
                d[i][j] = None

    for i in range(1, len(d)):
        d[i][0] = dividierteDifferenzenFormel(d, x, i, 0)
    return d

def dividierteDifferenzenFormel(d, x, n, m):
    if math.isnan(d[n][m+1]):
        d[n][m+1] = dividierteDifferenzenFormel(d, x, n, m+1)
    if math.isnan(d[n - 1][m]):
        d[n - 1][m] = dividierteDifferenzenFormel(d, x, n-1, m)
    return float(d[n][m + 1] - d[n - 1][m]) / (x[n] - x[m])


# Auswertung des Interpolationspolynoms nach dem Horner-aehnlichen Schema
def auswertung(x, d, x_werte):
    qk = 0.
    for i in range(0, len(d)):
        if i == 0:
            qk = d[len(d)-1][0]
        else:
            qk = d[len(d)-(i+1)][0] + (x - x_werte[len(d)-(i+1)]) * qk
    return qk


def getXs(m, f):
    x = np.empty(m, float)
    for i in range(0, m):
        x[i] = f(i, m)
    return x


def fillY(x, f):
    y = np.empty(len(x))
    for i in range(0, len(y)):
        y[i] = f(x[i])
    return y


def getSubplotPos(i):
    if i == 0:
        return 221
    elif i== 1:
        return 222
    elif i == 2:
        return 223


def zeroV(m):
    z = [0] * m
    return (z)


def cubic_spline(n, xn, a):
    h = zeroV(n - 1)
    alpha = zeroV(n - 1)

    l = zeroV(n + 1)
    u = zeroV(n)
    z = zeroV(n + 1)

    b = zeroV(n)
    c = zeroV(n + 1)
    d = zeroV(n)

    for i in range(n - 1):
        h[i] = xn[i + 1] - xn[i]

    for i in range(1, n - 1):
        alpha[i] = (3. / h[i]) * (a[i + 1] - a[i]) - (3. / h[i - 1]) * (a[i] - a[i - 1])

    l[0] = 1
    u[0] = 0
    z[0] = 0

    # II
    for i in range(1, n - 1):
        l[i] = 2 * (xn[i + 1] - xn[i - 1]) - h[i - 1] * u[i - 1]
        u[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    l[n] = 1
    z[n] = 0
    c[n] = 0

    for j in range(n - 2, -1, -1):
        c[j] = z[j] - u[j] * c[j + 1]
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3.
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    for j in range(n - 1):
        cub_graph(a[j], b[j], c[j], d[j], xn[j], xn[j + 1])

    plt.show()


def cub_graph(a, b, c, d, x_i, x_i_1):
    root = pylab.poly1d(x_i, True)
    poly = 0
    poly = d * (root) ** 3
    poly = poly + c * (root) ** 2
    poly = poly + b * root
    poly = poly + a

    pts = pylab.arange(x_i, x_i_1, 0.001)
    plt.plot(pts, poly(pts), '-')
    return


# Aufgabe 2a
# Beispiel aus dem Skript
x = [0, 1, 3]
y = [3, 2, 6]
d = dividierteDifferenzen (y, x)
print "d[0][0]:", d[0][0] # d_00 = 3
print "d[1][1]:", d[1][1] # d_11 = 2
print "d[2][2]:", d[2][2] # d_22 = 6
print "d[1][0]:", d[1][0] # d_10 = -1
print "d[2][1]:", d[2][1] # d_21 = 2
print "d[2][0]:", d[2][0] # d_20 = 1
# Auswertung an Stelle x = 2
print "Auswertung fuer x = 2:", auswertung(2, d, x)
# Auswertung an Stelle x = 3
print "Auswertung fuer x = 3:", auswertung(3, d, x)

#Aufgabe 2b Zur Visualisierung die Kommentarzeichen entfernen
f = function()
m = [7, 9, 11]

x = np.linspace(-5, 5, 100) # 100 linearly spaced numbers
y = np.empty(len(x), float)

for i in range(0, len(m)):
    x_interpol = getXs(m[i], f.mAufgabe2b)
    y_interpol = fillY(x_interpol, f.aufgabe2b)
    d = dividierteDifferenzen(y_interpol, x_interpol)
    for j in range(0, len(x)):
        y[j] = auswertung(x[j], d, x_interpol)
    #pylab.subplot(221+i)
    #pylab.plot(x, y)
    #pylab.title("m=" + str(m[i]))

y = fillY(x, f.aufgabe2b)

#pylab.subplot(224)
#pylab.plot(x, y)
#pylab.show()

#Aufgabe 2c Zur Visualisierung die Kommentarzeichen entfernen
# Die Interpolationsfunktionen sind von der originalen Funktion stark verschieden
# Dies aendert sich auch nicht mit zunehmenden ms
# Sollte die Grafikanzeige nicht funktionieren, ist eine PNG-Datei im Zip hinterlegt
x = np.linspace(-5, 5, 100) # 100 linearly spaced numbers
y = np.empty(len(x), float)

for i in range(0, len(m)):
    pos = getSubplotPos(i)
    x_interpol = getXs(m[i], f.mAufgabe2c)
    y_interpol = fillY(x_interpol, f.aufgabe2c)
    print x_interpol
    print y_interpol
    d = dividierteDifferenzen(y_interpol, x_interpol)
    for j in range(0, len(x)):
        y[j] = auswertung(x[j], d, x_interpol)
    #pylab.subplot(pos)
    #pylab.plot(x, y)
    #pylab.title("m=" + str(m[i]))

#pylab.subplot(224)
#pylab.plot(fillY(x, f.aufgabe2c), fillY(x, f.aufgabe2c))
#pylab.title("originale Funktion")
#pylab.show()


#Aufgabe 2d Sollte die Grafikanzeige nicht funktionieren, ist eine PNG-Datei im Zip hinterlegt
# Ab m = 8 hat man schon eine "gute" Naeherung erzielt
# Zur Visualisierung die Kommentarzeichen entfernen
m = [6, 7, 8]
t = np.linspace(0, 1, 100) # 100 linearly spaced numbers
x = np.empty(len(t), float)
y = np.empty(len(t), float)

for i in range(0, len(m)):
    pos = getSubplotPos(i)
    t_approx = getXs(m[i], f.mAufgabe2d) # t_approx sind die Stuetzstellen
    x_interpol = fillY(t_approx, f.aufgabe2dX)
    y_interpol = fillY(t_approx, f.aufgabe2dY)

    d_x = dividierteDifferenzen(x_interpol, t_approx)
    d_y = dividierteDifferenzen(y_interpol, t_approx)
    for j in range(0, len(t)):
        x[j] = auswertung(t[j], d_x, t_approx)
        y[j] = auswertung(t[j], d_y, t_approx)
    pylab.subplot(pos)
    pylab.plot(x, y)
    pylab.title("m="+str(m[i]))

pylab.subplot(224)
pylab.plot(fillY(t, f.aufgabe2dX), fillY(t, f.aufgabe2dY))
pylab.title("originale Funktion")
pylab.show()

#Aufgabe 4
#2b
m = [7, 9, 11]

#for i in range(0, len(m)):
    #x_interpol = getXs(m[i], f.mAufgabe2b)
    #y_interpol = fillY(x_interpol, f.aufgabe2b)
    #cubic_spline(m[i], x_interpol, y_interpol)
