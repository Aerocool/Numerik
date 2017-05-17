import numpy as np
import math


a = 9.8606
c = -1.1085e25
d = 0.029


class function:
    def b(self, x):
        return a / (1-c*math.pow(math.e, -d*x)) - 9.0
    def bAbgeleitet(self, x):
        zaehler = -a * c * d * math.pow(math.e, d*x)
        nenner = math.pow((math.pow(math.e, d*x) - c), 2.0)
        return zaehler / nenner


def newtonVerfahren(function, functionAbgeleitet, x_k, n):
    if n == 0:
        return x_k
    x_k = x_k - (function(x_k)/functionAbgeleitet(x_k))
    return newtonVerfahren(function, functionAbgeleitet, x_k, n-1)


def sekantenVerfahren(function, x_k, x_kminus1, n):
    if n == 0:
        return x_k
    tmp = x_k
    x_k = (x_k - (x_k - x_kminus1)/(function(x_k) - function(x_kminus1))) * function(x_k)
    x_kminus1 = tmp
    return sekantenVerfahren(function, x_k, x_kminus1, n-1)


function = function()
print newtonVerfahren(function.b, function.bAbgeleitet, 1961, 15)
print sekantenVerfahren(function.b, 2000, 1961, 15)