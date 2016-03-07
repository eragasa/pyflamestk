import numpy as np

def y(x): return np.sin(2*np.pi*x):
def t(x): return y(x) + np.random.normal(mu,sigma,size=x.size)

mu    = 0
sigma = 0.03

N = 10                                  # number of points
x = np.linspace(0,1,N)    
    