import os
import math

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
from fits import Fit, _simple_peak, linear_fit
from tka import TkaFile

from manuell import *
from number import *

os.chdir(os.path.dirname(os.path.realpath(__file__)))

SAVETYPE = "png"


def clean_filename(filename):
	return filename.replace(".", "_")


def local_fit(x, y, errors, mu, sigma):
	func = lambda x, mu, sigma, A, a0, a1: GAUSS(x, mu, sigma, A) + a0 + a1*x

	fit = Fit(func)
	fit.mu = _simple_peak(y, mu, sigma)
	fit.sigma = sigma
	fit.A = y[fit.mu]
	fit.a0 = 0
	fit.a1 = 0

	for i in range(5):
		lower = max(fit.mu - 5*fit.sigma, 0)
		upper = min(fit.mu + 5*fit.sigma, len(x))
		lx = x[lower:upper]
		ly = y[lower:upper]
		lerrors = errors[lower:upper]
		if i <= 2:
			fit.fit(lx, ly, lerrors, solver="curve_fit")
		else:
			fit.fit(lx, ly, lerrors, solver="minuit", method="migrad")

	return fit

def plot_energy_spectra():
	mus = []
	widths = []
	energies = []

	total_count = 0
	total_time = 0

	for filename, meta in energy_spectra.items():
		tka = TkaFile("data/" + filename)
		x = np.arange(len(tka))
		y = tka.data

		errors = np.sqrt(y+1)

		plt.clf()
		plt.plot(x, y, '.', color="black")

		for i, (mu0, sigma0, energy, sigma_energy) in enumerate(meta["peaks"], 1):
			fit = local_fit(x,y,errors, mu=mu0, sigma=sigma0)
			fit.plot(fit.mu-5*fit.sigma, fit.mu+5*fit.sigma, 1000, zorder=10000, color="red")
			mus.append((fit.mu, fit.error_mu))
			widths.append((fit.sigma, fit.error_sigma))
			energies.append((energy, sigma_energy))
			print(meta["element"], i, fit)

		plt.xlabel("Kanal")
		plt.ylabel("Count")
		plt.xlim(0, 2**14)
		plt.title(meta["title"])
		plt.savefig("out/" + clean_filename(filename) + "_all." + SAVETYPE)

		total_count += y.sum()
		total_time += tka.real_time

	mus, error_mus = np.array(mus).T
	widths, error_widths = np.array(widths).T
	energies, error_energies = np.array(energies).T

	# Kalibration
	plt.clf()
	fit = linear_fit(mus, energies, xerr=error_mus, yerr=error_energies, slope=0.1, offset=-20)
	plt.errorbar(mus, energies, xerr=error_mus, yerr=error_energies, fmt=',')
	fit.plot(mus.min(), mus.max(), box='tl', units={"slope": "eV / Channel", "offset": "eV"}, factors={"slope": 1000, "offset": 1000})
	plt.xlabel("Kanal")
	plt.ylabel("Energie / keV")
	plt.savefig("out/calib_fit." + SAVETYPE)

	plt.clf()
	errs = fit.combine_errors(mus, xerr=error_mus, yerr=error_energies)
	fit.plot_residual(mus, energies, errs, box='tl', fmt=',')
	plt.xlabel("Kanal")
	plt.ylabel("Energie / keV")
	plt.savefig("out/calibresiduum." + SAVETYPE)


	# Energieauflösung
	X = energies
	Y = np.power(widths, 2)
	SX = error_energies
	SY = 2*widths*error_widths

	func = lambda x, a0, a1,: a0+x*a1

	fit = Fit(func)
	fit.a0 = 1
	fit.a1 = 1

	for _ in range(10):
		err = fit.combine_errors(X, SX, SY)
		fit.fit(X, Y, err)

	plt.clf()
	plt.errorbar(X, Y, xerr=SX, yerr=SY, fmt=',')
	fit.plot(X.min(), X.max(), box='br', units={'a0': 'keV', 'a1': r'\sqrt{keV}'})
	plt.xlabel(r"$ E $ / keV")
	plt.ylabel(r"$ \Delta E^2 / keV^2$")
	plt.title(r'Energieauflösung: Fit zu $ \Delta E^2 = a_0 + a_1 E $')
	plt.savefig("out/energyresolution_fit." + SAVETYPE)

	if fit.a1 >= 0:
		a = math.sqrt(fit.a1)
		error_a = fit.error_a1 / (2*a)
		print("a =", formatQuantity(a, error_a))
	else:
		print("Parameter a complex.")

	if fit.a0 >= 0:
		b = math.sqrt(fit.a0)
		error_b = fit.error_a0 / (2*b)
		print("b =", formatQuantity(b, error_b))
	else:
		print("Parameter b complex.")



	plt.clf()
	err = fit.combine_errors(X, SX, SY)
	fit.plot_residual(X, Y, err, box='tr', fmt=",")
	plt.xlabel(r"$ E $ / keV")
	plt.ylabel(r"$ \Delta E^2 / keV^2$")
	plt.title("Energieauflösung: Residuen")
	plt.savefig("out/energyresolution_residuum." + SAVETYPE)

	print("Total number of events:", total_count)
	print("Total measuring time:", total_time, "s")
	print("Event rate:", (total_count/total_time), "1/s")


def plot_sca_spectrum(filename, stepsize):
	y = np.loadtxt(filename)
	y[0] = 0
	x = stepsize * np.arange(len(y))
	plt.plot(x,y, 's')
	plt.xlabel("Window Offset / V")
	plt.ylabel("Count")
	plt.show()

def plot_delay_spectrum():
	x1, y1 = np.loadtxt("data/pm1_delay.txt", unpack=True)
	x2, y2 = np.loadtxt("data/pm2_delay.txt", unpack=True)

	time = np.hstack((-x1, x2))
	count = np.hstack((y1, y2))
	error_count = np.sqrt(count + 1)
	plt.errorbar(time, count, yerr=error_count, fmt='s')

	fit = Fit(GAUSS)
	fit.mu = -45
	fit.sigma = 60
	fit.fit(time, count, error_count)
	fit.plot(time.min(), time.max(), box="tr", units={"mu": "ns", "sigma": "ns"})
	print(fit)

	plt.xlabel(r"$\Delta T / \textrm{ns}$")
	plt.ylabel("Count")
	plt.show()

def plot_coincidence_na():
	angle, count1, count2 = np.loadtxt("data/natrium_winkel.txt")


if __name__=="__main__":
	#plot_energy_spectra()
	#plot_delay_spectrum()
	#plot_sca_spectrum("data/pm1_sca_windowsearch2.txt", 0.2)
	#plot_sca_spectrum("data/pm2_sca_windowsearch2.txt", 0.2)
	plot_coincidence_na()
