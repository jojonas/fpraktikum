import math
import numpy as np
from scipy.optimize import curve_fit

from funcs import GAUSS

def chisqndf(x, y, fit, yerr, nparams=0):
	return (np.power((fit(x) - y)/yerr, 2)).sum() / (len(x)-nparams)

def easy_peak(data, index):
	for i in range(len(data)):
		if index+1 < len(data) and data[index+1] > data[index]:
			index += 1
		elif index-1 >= 0 and data[index-1] > data[index]:
			index -= 1
		else:
			break
	return index

def gauss_fit(data, mu=0, sigma=10, A=100, N=5, errors=None):
	peak = easy_peak(data, mu)
	for i in range(N):
		lower = max(int(peak-sigma), 0)
		upper = min(int(peak+sigma), len(data))

		xdata = np.arange(lower, upper)
		ydata = data[lower:upper]
		yerrors = errors[lower:upper]
		popt, pcov = curve_fit(GAUSS, xdata, ydata, p0=(A,mu,sigma), sigma=yerrors)
		A, mu, sigma = popt
		sigma = np.abs(sigma)

	return mu, sigma, A
