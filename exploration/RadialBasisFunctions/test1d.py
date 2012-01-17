import numpy as num
import matplotlib
import matplotlib.pyplot as plt
import math
from scipy.integrate import trapz

def kernel(sigma):
    sz = 6*num.ceil(sigma)+1
    x = num.linspace(-3*sigma, 3*sigma, sz);
    return matplotlib.mlab.normpdf(x, 0, sigma)
    
def makeI1():
    phantom = num.zeros(1000)
    phantom[250:750] = 1;
    return num.convolve(phantom, kernel(10), 'same')

def makeI0(offset):
    phantom = num.zeros(1000)
    phantom[(250+offset):(750+offset)] = 1;
    return num.convolve(phantom, kernel(10), 'same')

def csrbf(n, l, r):
    """compactly supported radial basis function of degree one"""
    p = num.floor(n/2) + l + 1
    return num.power((num.maximum(0, 1-r)), p+1)*((p+1)*r+1)

def phi(x, alpha, beta, gamma):
    """paramterized levelset function"""
    def f(xp, a, b, g):
        delta = xp - g
        return a*csrbf(1,1,abs(b*delta))

    vphi = num.vectorize(f);
    n = num.size(alpha)
    sum = num.zeros(num.size(x));
    for j in range(0,n):
        sum += vphi(x, alpha[j], beta[j], gamma[j])
    return sum

def heaviside(x, eps):
    """A continous approximation to the Heaviside function"""
    def f(xp):
        return 0.5*(1 + (2/math.pi)*math.atan(math.pi*xp/eps))

    vh = num.vectorize(f)
    return vh(x)

def f(I0, I1, tau, alpha, beta, gamma, lam, mu):
    """The objective function"""
    (x, step) = num.linspace(-1, 1, 1000, retstep=True)
    def integrand1(xp):
        return tau*tau*heaviside(phi(xp,alpha, beta,gamma), 0.01)

    def integrand2(xp):
        moving = num.interp(xp-tau, xp, I0)
        residual = I1 - I0;
        return heaviside(phi(xp,alpha, beta,gamma),0.01)*num.power(residual, 2)
        
    def integrand3(xp):
        return heaviside(phi(xp,alpha, beta,gamma), 0.01)
        
    return trapz(integrand1(x)+lam*integrand2(x)+mu*integrand3(x), dx=step)

                                 
#    term1 = num.add.reduce(tau*tau*vh(phi(x, alpha, beta, gamma),0.01))
#    moving = interpolate(x-tau,I0)
#    residual = num.power(I1 - moving,2)
#    term2 = num.add.reduce(vh(phi(x, alpha, beta, gamma),0.01))
    
def test():
    x = num.linspace(-1, 1, 1000)

    I1 = makeI1()
    I0 = makeI0(10)

    a = num.zeros(1)
    B = 4*num.ones(1)
    g = -0.5*num.ones(1)
    
    tau = 0
    lam = 1
    mu = 1

    print f(I0, I1, tau, a, B, g, lam, mu)

    plt.plot(x, phi(x, a, B, g))
    plt.show()
    
    
