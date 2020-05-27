
import scipy 

import numpy as np

a = 0

b = 2


def f1(x):
    return x**2-2


print (scipy.optimize.brentq(f1, a, b))

def f2(x):
    return x**5 -6.6*x**4+5.12*x**3+21.312*x**2-38.016*x+17.28

print (scipy.optimize.brentq(f2, a, b))

def f3(x):
    return (x-1.5)* np.exp(-4*(x-1.5)**2)

print (scipy.optimize.brentq(f3, a, b))