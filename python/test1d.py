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
    return num.power((num.maximum(0, 1-r)), 3)*(3*r+1)

### @export "dcsrbf"
def dcsrbf(r):
    """derivative compactly supported radial basis function of degree one"""
    return 3*(num.power(num.maximum(0, 1-r), 3) - num.power(num.maximum(0, 1-r),2)*(3*r+1))

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
    def f(xp):
        delta = xp - gamma_j
        return csrbf(math.sqrt(beta_j*delta*beta_j*delta))
    vf = num.vectorize(f)
    return vf(x)

def dphidbeta_j(x, alpha_j, beta_j, gamma_j):
    def f(xp):
        delta = xp - gamma_j
        r = math.sqrt(beta_j*delta*beta_j*delta) + 0.001
        return alpha_j*dcsrbf(r)*(beta_j*delta*delta)/r
    vf = num.vectorize(f)
    return vf(x)

def dphidgamma_j(x, alpha_j, beta_j, gamma_j):
    def f(xp):
        delta = xp - gamma_j
        r = math.sqrt(beta_j*delta*beta_j*delta) + 0.001
        return alpha_j*dcsrbf(r)*(-beta_j*beta_j*delta)/r
    vf = num.vectorize(f)
    return vf(x)

### @export "H"
def heaviside(x, eps=0.001):
    """A continuous approximation to the Heaviside function"""
    def f(xp):
        return 0.5*(1 + (2/math.pi)*math.atan(math.pi*xp/eps))
    vh = num.vectorize(f)
    return vh(x)

def deltaH(x, eps=0.001):
    """A continuous approximation to the Heaviside function"""
    def f(xp):
        return (eps/(eps*eps + math.pi*math.pi*xp*xp))
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
    def integrand2(xp):
        moving = num.interp(xp-tau, xp, I0)
        residual = I1 - moving;
        return heaviside(-phi(xp,alpha, beta,gamma))*num.power(residual, 2)
    def integrand3(xp):
        return heaviside(phi(xp,alpha, beta,gamma))
    return trapz(lam*integrand2(x)+mu*integrand3(x), dx=step)
    
def dfdtau(I0, I1, tau, alpha, beta, gamma, lam, mu):
    (x, step) = num.linspace(-1, 1, 1000, retstep=True)
    def integrand2(xp):
        moving = num.interp(xp-tau, xp, I0)
        residual = I1 - moving;
        return 2*tau*heaviside(-phi(xp,alpha, beta,gamma))*residual*num.gradient(moving, step)
    return trapz(lam*integrand2(x), dx=step)

def dfdalpha_j(I0, I1, tau, alpha, beta, gamma, lam, mu, j):
    (x, step) = num.linspace(-1, 1, 1000, retstep=True)
    def integrand2(xp):
        moving = num.interp(xp-tau, xp, I0)
        residual = I1 - moving;
        return -num.power(residual, 2)*deltaH(phi(xp,alpha, beta,gamma))*dphidalpha_j(xp, alpha[j], beta[j], gamma[j])
    def integrand3(xp):
        return deltaH(phi(xp,alpha, beta,gamma))*dphidalpha_j(xp, alpha[j], beta[j], gamma[j])
    return trapz(lam*integrand2(x)+mu*integrand3(x), dx=step)

def dfdbeta_j(I0, I1, tau, alpha, beta, gamma, lam, mu, j):
    (x, step) = num.linspace(-1, 1, 1000, retstep=True)
    def integrand2(xp):
        moving = num.interp(xp-tau, xp, I0)
        residual = I1 - moving;
        return -num.power(residual, 2)*deltaH(phi(xp,alpha, beta,gamma))*dphidbeta_j(xp, alpha[j], beta[j], gamma[j])
    def integrand3(xp):
        return deltaH(phi(xp,alpha, beta,gamma))*dphidbeta_j(xp, alpha[j], beta[j], gamma[j])
    return trapz(lam*integrand2(x)+mu*integrand3(x), dx=step)

def dfdgamma_j(I0, I1, tau, alpha, beta, gamma, lam, mu, j):
    (x, step) = num.linspace(-1, 1, 1000, retstep=True)
    def integrand2(xp):
        moving = num.interp(xp-tau, xp, I0)
        residual = I1 - moving;
        return -num.power(residual, 2)*deltaH(phi(xp,alpha, beta,gamma))*dphidgamma_j(xp, alpha[j], beta[j], gamma[j])
    def integrand3(xp):
        return deltaH(phi(xp,alpha, beta,gamma))*dphidgamma_j(xp, alpha[j], beta[j], gamma[j])
    return trapz(lam*integrand2(x)+mu*integrand3(x), dx=step)

### @export "wrapper"
def wrapper(theta, I0, I1, lam, mu):
    return f(I0, I1, theta[0], theta[1:8], theta[8:15], theta[15:23], lam, mu)

def gradf(theta, I0, I1, lam, mu):
    tau=theta[0]
    alpha = theta[1:8]
    beta = theta[8:15]
    gamma = theta[15:23]
    grad = num.zeros(num.size(theta))
    grad[0] = dfdtau(I0, I1, tau, alpha, beta, gamma, lam, mu)
    for j in range(0, num.size(alpha)):
        grad[j+1] = dfdalpha_j(I0, I1, tau, alpha, beta, gamma, lam, mu, j)
        grad[j+8] = dfdbeta_j(I0, I1, tau, alpha, beta, gamma, lam, mu, j)
        grad[j+15] = dfdgamma_j(I0, I1, tau, alpha, beta, gamma, lam, mu, j)
    return grad

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
    theta0[0] = 0 #initial translation
    theta0[1] = -0.2 #initial alpha_1
    theta0[2] = 0.2 #initial alpha_2
    theta0[3] = -0.2 #initial alpha_3
    theta0[4] = 0.2 #initial alpha_4
    theta0[5] = -0.2 #initial alpha_5
    theta0[6] = 0.2 #initial alpha_6
    theta0[7] = -0.2 #initial alpha_7
                
    theta0[8:15] = 2.0 #initial beta_j

    theta0[15:23] = num.linspace(-1,1, m0) #initial gamma_j

    print wrapper(theta0, M, T, lam, mu)
    print gradf(theta0, M, T, lam, mu)
    
#    df = gradf(theta0, M, T, lam, mu)
#    print num.size(theta0), " ", num.size(df)
#    print df
    
#    opt = theta0
    # opt = fmin_bfgs(wrapper,
    #                 theta0,
    #                 gradf,
    #                 args=(M, T, lam, mu),
    #                 callback=monitor,
    #                 maxiter=50)

    # opt = fmin_bfgs(wrapper,
    #                 theta0,
    #                 args=(M, T, lam, mu),
    #                 callback=monitor,
    #                 maxiter=20)

    # print "start"
    # print "     params = ", theta0
    # print "     grad = ", gradf(theta0, M, T, lam, mu)
    # print "     objective = ", wrapper(theta0, M, T, lam, mu)
    # print "     ROI size = ", trapz(heaviside(phi(x,theta0[1:8], theta0[8:15], theta0[15:23])), dx=step)
    # print "end"
    # print "   params = ", opt
    # print "   grad = ", gradf(opt, M, T, lam, mu)
    # print "   objective = ", wrapper(opt, M, T, lam, mu)
    # print "   ROI size = ", trapz(heaviside(phi(x,opt[1:8], opt[8:15], opt[15:23])), dx=step)
    
    # plt.figure(1)
    # plt.plot(x, heaviside(phi(x, theta0[1:8], theta0[8:15], theta0[15:23])))
    # plt.xlabel('x')

    # plt.figure(2)
    # plt.plot(x, heaviside(phi(x, opt[1:8], opt[8:15], opt[15:23])))
    # plt.xlabel('x')



    # plt.show()


    
    
