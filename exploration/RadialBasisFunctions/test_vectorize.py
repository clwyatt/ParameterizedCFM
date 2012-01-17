import numpy as num
import math
import matplotlib.pyplot as plt

def heaviside(x, eps):

    def f(xp):
        return 0.5*(1 + (2/math.pi)*math.atan(math.pi*xp/eps))

    vh = num.vectorize(f)
    return vh(x)
    
x = num.linspace(-1, 1, 1000)
y = heaviside(x, 0.01)

plt.plot(x,y)
plt.show()



