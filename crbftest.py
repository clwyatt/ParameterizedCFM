import math
import numpy as num
import matplotlib.pyplot as plt

def csrbf(r):
    """compactly supported radial basis function of degree one"""
    def f(rp):
        return num.power((num.maximum(0, 1-rp)), 3)*(3*rp+1)
    fvec = num.vectorize(f)
    return fvec(r)

def dcsrbf(r):
    """derivative compactly supported radial basis function of degree one"""
    def f(rp):
        return 3*(num.power(num.maximum(0, 1-rp), 3) - num.power(num.maximum(0, 1-rp),2)*(3*rp+1))
    fvec = num.vectorize(f)
    return fvec(r)

def test():
    (r, step) = num.linspace(0, 2, 1000, retstep=True)

    plt.plot(r, csrbf(r))
    plt.plot(r, dcsrbf(r))
    plt.xlabel('r')
    plt.show()

