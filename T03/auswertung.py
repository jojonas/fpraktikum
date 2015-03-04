import os

from pylab import rcParams
rcParams['savefig.dpi'] = 100
rcParams['text.usetex'] = True
rcParams['text.latex.unicode'] = True
#rcParams['text.latex.preamble']=[r'\usepackage{lmodern}']
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'Computer Modern'
rcParams['font.size'] = 18
rcParams['savefig.bbox'] = 'tight'
rcParams['axes.grid'] = True

import numpy as np
import matplotlib.pyplot as plt

from funcs import *
from fits import Fit
from tka import TkaFile

os.chdir(os.path.dirname(os.path.realpath(__file__)))

def local_fit(x, y, errors, mu, sigma, A):
	func = lambda x, mu, sigma, A, a0, a1: GAUSS(x, mu, sigma, A) + a0 + a1*x

	fit = Fit(func)
	fit.mu = mu
	fit.sigma = sigma
	fit.A = A
	fit.a0 = 0
	fit.a1 = 0

	lower = max(fit.mu - 5*fit.sigma, 0)
	upper = min(fit.mu + 5*fit.sigma, len(x))
	lx = x[lower:upper]
	ly = y[lower:upper]
	lerrors = errors[lower:upper]
	fit.fit(lx, ly, lerrors)

	return fit

def show(filename):
	tka = TkaFile(filename)
	x = np.arange(len(tka))
	y = tka.data

	print("Frequenz:", (y.sum() / tka.live_time))

	errors = np.sqrt(y)
	errors[errors==0] = 1

	# fit = local_fit(x,y,errors, mu=5400, sigma=200, A=3000)
	# fit.plot(fit.mu-5*fit.sigma, fit.mu+5*fit.sigma, 1000, zorder=10000)
	# print(fit.mu)
	#
	# fit = local_fit(x,y,errors, mu=13050, sigma=400, A=500)
	# fit.plot(fit.mu-5*fit.sigma, fit.mu+5*fit.sigma, 1000, zorder=10000)
	# print(fit.mu)

	plt.plot(x, y, '.')
	plt.xlabel("Kanal")
	plt.ylabel("Count")
	plt.xlim(0, 2**14)
	plt.show()

def plot_sca_spectrum(filename):
	y = np.loadtxt(filename)
	y[0] = 0
	x = np.arange(len(y)) / 10.0
	plt.plot(x,y, 's')
	plt.show()

def plot_delay_spectrum():
	x1, y1 = np.loadtxt("data/pm1_delay.txt", unpack=True)
	x2, y2 = np.loadtxt("data/pm2_delay.txt", unpack=True)
	x1 = -x1
	plt.errorbar(x1,y1, yerr=np.sqrt(y1), fmt='s')
	plt.errorbar(x2,y2, yerr=np.sqrt(y2), fmt='s')
	plt.xlabel("$\Delta T$")
	plt.ylabel("Count")
	plt.show()

if __name__=="__main__":
	#show("data/EnergiespektrumNa.TKA")
	#plot_sca_spectrum("data/pm1_sca_windowsearch.txt")
	#plot_sca_spectrum("data/pm2_sca_windowsearch.txt")
	plot_delay_spectrum()
