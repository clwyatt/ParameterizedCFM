### @export "imports"
import math
import numpy as num
import matplotlib.pyplot as plt

### @export "csrbf"
def csrbf(r):
    """compactly supported radial basis function of degree one"""
    return num.power((num.maximum(0, 1-r)), 3)*(3*r+1)

### @export "levelset"
def phi(x, alpha, beta, gamma):
    """paramterized levelset function"""
    def f(xp, a, b, g):
        delta = xp - g
        return a*csrbf(math.sqrt(b*delta*b*delta + 0.001))
    vphi = num.vectorize(f);
    n = num.size(alpha)
    sum = num.zeros(num.size(x));
    for j in range(0,n):
        sum += vphi(x, alpha[j], beta[j], gamma[j])
    return sum

### @export "H"
def heaviside(x, eps=0.001):
    """A continuous approximation to the Heaviside function"""
    def f(xp):
        return 0.5*(1 + (2/math.pi)*math.atan(math.pi*xp/eps))
    vh = num.vectorize(f)
    return vh(x)


if __name__ == "__main__":
    opt = num.genfromtxt('result.dat')
    (x, step) = num.linspace(-1, 1, 1000, retstep=True)
    plt.plot(x, heaviside(phi(x, opt[1:8], opt[8:15], opt[15:23])))
    plt.xlabel('x')
    
    plt.show()


    
    
