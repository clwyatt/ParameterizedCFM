import numpy as num
from scipy.integrate import trapz

def afunction(x, param):
    return x + param

def test():
    (x, step) = num.linspace(-1,1,1000, retstep=True)
    y = afunction(x, 1);
    print trapz(y,dx=step)
    
