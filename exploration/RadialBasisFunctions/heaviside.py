import numpy as num
import matplotlib.pyplot as plt
import math

def heaviside(x, eps):
    return 0.5*(1 + (2/math.pi)*math.atan(math.pi*x/eps))
    
x = num.linspace(-1, 1, 1000)

vh = num.vectorize(heaviside)

plt.plot(vh(x,0.01))
plt.show()
