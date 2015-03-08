import os

from pylab import rcParams
rcParams['savefig.dpi'] = 100
#rcParams['text.usetex'] = True
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
from fits import Fit, _simple_peak
from tka import TkaFile

from manuell import *

os.chdir(os.path.dirname(os.path.realpath(__file__)))

def local_fit(x, y, errors, mu, sigma):
	func = lambda x, mu, sigma, A, a0, a1: GAUSS(x, mu, sigma, A) + a0 + a1*x

	fit = Fit(func)
	fit.mu = _simple_peak(y, mu, sigma)
	fit.sigma = sigma
	fit.A = y[fit.mu]
	fit.a0 = 0
	fit.a1 = 0

	lower = max(fit.mu - 5*fit.sigma, 0)
	upper = min(fit.mu + 5*fit.sigma, len(x))
	lx = x[lower:upper]
	ly = y[lower:upper]
	lerrors = errors[lower:upper]
	fit.fit(lx, ly, lerrors)

	return fit

def plot_energy_spectra():
	for filename, meta in energy_spectra.items():
		tka = TkaFile("data/" + filename)
		x = np.arange(len(tka))
		y = tka.data

		errors = np.sqrt(y+1)
		plt.plot(x, y, '.', color="black")

		for i, (mu0, sigma0) in enumerate(meta["peaks"], 1):
			fit = local_fit(x,y,errors, mu=mu0, sigma=sigma0)
			fit.plot(fit.mu-5*fit.sigma, fit.mu+5*fit.sigma, 1000, zorder=10000, color="red")
			print(meta["element"], i, fit)

		plt.xlabel("Kanal")
		plt.ylabel("Count")
		plt.xlim(0, 2**14)
		plt.title(meta["title"])
		plt.show()

def plot_sca_spectrum(filename, stepsize):
	y = np.loadtxt(filename)
	y[0] = 0
	x = stepsize * np.arange(len(y))
	plt.plot(x,y, 's')
	plt.show()

def plot_delay_spectrum():
	x1, y1 = np.loadtxt("data/pm1_delay.txt", unpack=True)
	x2, y2 = np.loadtxt("data/pm2_delay.txt", unpack=True)
	x1 = -x1
	plt.errorbar(x1,y1, yerr=np.sqrt(y1), fmt='s')
	plt.errorbar(x2,y2, yerr=np.sqrt(y2), fmt='s')
	plt.xlabel(r"$\Delta T / \textrm{ns}$")
	plt.ylabel("Count")
	plt.show()

if __name__=="__main__":
	plot_energy_spectra()
	#plot_sca_spectrum("data/pm1_sca_windowsearch2.txt", 0.2)
	#plot_sca_spectrum("data/pm2_sca_windowsearch2.txt", 0.2)
	#plot_delay_spectrum()
