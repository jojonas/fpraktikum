import numpy as np

GAUSS = lambda x, A, mu, sigma: A*np.exp(-np.power(x-mu, 2)/(2*sigma**2))
