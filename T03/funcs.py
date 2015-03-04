import numpy as np

GAUSS = lambda x, mu=1, sigma=1, A=1: A*np.exp(-np.power(x-mu, 2)/(2*sigma**2))
LINEAR = lambda x, slope=1, offset=0: slope*x + offset
SQRT = lambda x, A=1, offset=0: A*np.sqrt(x)+offset
