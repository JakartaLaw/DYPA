# Numerical intergration
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.hermite import hermgauss


# %%
# Monte carlo integration

def f():
    """function to be approximated"""

    return np.random.normal()**2

# %%
def monte_carlo_integration(iterations, func):
    """Monte Carlo Integration"""
    return np.mean([func() for i in range(iterations)])

# %%
avg = monte_carlo_integration(100000, f)
print(avg)

# %%
# Equiprobale integration

class struct():
    pass

def Equiprobale(S, func, inv_func):
    """
    S: number of points in grid,
    func : function for numerical integration (pdf)
    inv_func : function inversed
    """

    grid = np.linspace(0, 1, S)
    z = map(inv_func, grid)

# %%
# Gauss Hermite

def _f(x):
    return x**2

mu, sigma = 0, 1
N_GH = 9

_x, _w = hermgauss(N_GH)

x = _x*np.sqrt(2) * sigma + mu
w = _w/np.sqrt(np.pi)

np.sum(_f(x)*w)
