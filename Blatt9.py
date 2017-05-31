import numpy as np
import math
import matplotlib.pyplot as plt
import pylab

class function:
    def aufgabe1(self, x):
        return float((-16.0/15.0) + (69.0)/(5.0)*x**2)


def fillY(x, f):
    y = np.empty(len(x))
    for i in range(0, len(y)):
        y[i] = f(x[i])
    return y

x_i = [-2, 0, 1, 2]
y_i = [4, 0, 1, -4]

x = np.linspace(-2, 2, 100) # 100 linearly spaced numbers
y = fillY(x, function().aufgabe1)

pylab.plot(x, y)
pylab.plot(x_i, y_i, 'ro')
pylab.title("Aufgabe 1")
pylab.show()