import numpy as np

GAUSS = lambda x, mu=1, sigma=1, A=1: A*np.exp(-np.power(x-mu, 2)/(2*sigma**2))
LINEAR = lambda x, slope=1, offset=0: slope*x + offset
SQRT = lambda x, A=1, offset=0: A*np.sqrt(x)+offset
CONSTANT = lambda x, c=0: c
EXPONENTIAL_DECAY = lambda x, A, T, offset: A*np.exp(-x/T)+offset
LORENTZ = CAUCHY = lambda x, s, t: s / np.pi / (np.power(s,2) + np.power(x-t, 2))

POLY0 = lambda x, a0=0: a0*np.ones_like(x)
POLY1 = lambda x, a0=0, a1=1: a0 + x*a1
POLY2 = lambda x, a0=0, a1=1, a2=1: a0 + x*a1 + np.power(x, 2)*a2
POLY3 = lambda x, a0=0, a1=1, a2=1, a3=1: a0 + x*a1 + np.power(x, 2)*a2 + np.power(x, 3)*a3
POLY4 = lambda x, a0=0, a1=1, a2=1, a3=1, a4=1: a0 + x*a1 + np.power(x, 2)*a2 + np.power(x, 3)*a3 + np.power(x, 4)*a4
POLY5 = lambda x, a0=0, a1=1, a2=1, a3=1, a4=1, a5=1: a0 + x*a1 + np.power(x, 2)*a2 + np.power(x, 3)*a3 + np.power(x, 4)*a4 + np.power(x, 5)*a5
POLY6 = lambda x, a0=0, a1=1, a2=1, a3=1, a4=1, a5=1, a6=1: a0 + x*a1 + np.power(x, 2)*a2 + np.power(x, 3)*a3 + np.power(x, 4)*a4 + np.power(x, 5)*a5 + np.power(x, 6)*a6
