### @export "imports"
import math
import numpy as num
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf 
from scipy.integrate import trapz
from scipy.optimize import fmin_bfgs

### @export "kernel"
def kernel(sigma):
    sz = 6*num.ceil(sigma)+1
    x = num.linspace(-3*sigma, 3*sigma, sz);
    return normpdf(x, 0, sigma)

### @export "phantom1"
def makeI1():
    phantom = num.zeros(1000)
    phantom[250:750] = 1;
    return num.convolve(phantom, kernel(10), 'same')

### @export "phantom2"
def makeI2():
    phantom = num.zeros(1000)
    phantom[250:750] = 1
    phantom[490:510] = 1.2
    return num.convolve(phantom, kernel(10), 'same')

### @export "phantom0"
def makeI0(offset):
    phantom = num.zeros(1000)
    phantom[(250+offset):(750+offset)] = 1;
    return num.convolve(phantom, kernel(10), 'same')

### @export "csrbf"
def csrbf(r):
    """compactly supported radial basis function of degree one"""
    return num.power((num.maximum(0, 1-rp)), 2)*(2*rp+1)

### @export "dcsrbf"
def dcsrbf(r):
    """derivative compactly supported radial basis function of degree one"""
    return -6*rp*(num.maximum(0, 1-rp))

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

def dphidalpha_j(x, alpha_j, beta_j, gamma_j):
    """derivative of phi wrt alpha_j"""
    delta = x - gamma_j
    return csrbf(math.sqrt(beta_j*delta*beta_j*delta))

def dphidbeta_j(x, alpha_j, beta_j, gamma_j):
    delta = x - gamma_j
    r = math.sqrt(beta_j*delta*beta_j*delta) + 0.001
    return alpha_j*dcsrbf(r)*(beta_j*delta*delta)/r

def dphidgamma_j(x, alpha_j, beta_j, gamma_j):
    delta = x - gamma_j
    r = math.sqrt(beta_j*delta*beta_j*delta) + 0.001
    return alpha_j*dcsrbf(r)*(-beta_j*beta_j*delta)/r

### @export "H"
def heaviside(x, eps=0.001):
    """A continuous approximation to the Heaviside function"""
    def f(xp):
        return 0.5*(1 + (2/math.pi)*math.atan(math.pi*xp/eps))
    vh = num.vectorize(f)
    return vh(x)

#def heaviside(x, eps=0.001):
#    """A continuous approximation to the Heaviside function"""
#    def f(xp):
#        if(xp > eps): return 1
#        if(xp < eps): return 0
#        return 0.5*(1 + xp/eps + (1/math.pi)*math.sin(math.pi*xp/eps))
#    vh = num.vectorize(f)
#    return vh(x)

### @export "objective"
def f(I0, I1, tau, alpha, beta, gamma, lam, mu):
    """The objective function"""
    (x, step) = num.linspace(-1, 1, 1000, retstep=True)
    def integrand1(xp):
        return tau*tau*heaviside(-phi(xp,alpha, beta,gamma))
    def integrand2(xp):
        moving = num.interp(xp-tau, xp, I0)
        residual = I1 - moving;
        return heaviside(-phi(xp,alpha, beta,gamma))*num.power(residual, 2)
    def integrand3(xp):
        return heaviside(phi(xp,alpha, beta,gamma))
    return trapz(lam*integrand2(x)+mu*integrand3(x), dx=step)

### @export "wrapper"
def wrapper(theta, I0, I1, lam, mu):
    return f(I0, I1, theta[0], theta[1:8], theta[8:15], theta[15:23], lam, mu)

def monitor(x):
    """Callback during optimization"""
    print x

### @export "test"
def test():
    (x, step) = num.linspace(-1, 1, 1000, retstep=True)

    T = makeI2()
    M = makeI0(100)

    lam = 10
    mu = 0.01

    m0 = 7
    theta0 = num.zeros(1+m0*3)
    theta0[0] = -step*100 #initial translation
    theta0[1] = -0.2 #initial alpha_1
    theta0[2] = 0.2 #initial alpha_2
    theta0[3] = -0.2 #initial alpha_3
    theta0[4] = 0.2 #initial alpha_4
    theta0[5] = -0.2 #initial alpha_5
    theta0[6] = 0.2 #initial alpha_6
    theta0[7] = -0.2 #initial alpha_7
                
    theta0[8:15] = 2.0 #initial beta_j

    theta0[15:23] = num.linspace(-1,1, m0) #initial gamma_j

#    opt = theta0
    opt = fmin_bfgs(wrapper,
                    theta0, args=(M, T, lam, mu),
                    callback=monitor,
                    maxiter=20)

    print "start"
    print "     params = ", theta0
    print "     objective = ", wrapper(theta0, M, T, lam, mu)
    print "     ROI size = ", trapz(heaviside(phi(x,theta0[1:8], theta0[8:15], theta0[15:23])), dx=step)
    print "end"
    print "   params = ", opt
    print "   objective = ", wrapper(opt, M, T, lam, mu)
    print "   ROI size = ", trapz(heaviside(phi(x,opt[1:8], opt[8:15], opt[15:23])), dx=step)
    
    plt.figure(1)
    plt.plot(x, heaviside(phi(x, theta0[1:8], theta0[8:15], theta0[15:23])))
    plt.xlabel('x')

    plt.figure(2)
    plt.plot(x, heaviside(phi(x, opt[1:8], opt[8:15], opt[15:23])))
    plt.xlabel('x')



    plt.show()


    
    
