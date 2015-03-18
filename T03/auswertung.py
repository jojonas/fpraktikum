import os
import math

from pylab import rcParams
rcParams['savefig.dpi'] = 200
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
from uncertainties import ufloat

from funcs import *
from fits import Fit, _simple_peak, linear_fit
from tka import TkaFile
from latextable import LatexTable

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


	with LatexTable("out/calib_peaks.tex") as peaktable, LatexTable("out/rates.tex") as ratetable:
		peaktable.header("El.", "Höhe", r"Channel $x_\textrm{max}$", r"Breite $\Delta x$", r'$\chi^2/\textrm{ndf}$', "Energie / keV", lineafter=1)
		ratetable.header("El.", "Events", "Messzeit", "Rate")

		for filename, meta in energy_spectra.items():
			peaktable.hline()

			tka = TkaFile("data/" + filename)
			x = np.arange(len(tka))
			y = tka.data

			errors = np.sqrt(y+1)

			plt.clf()
			plt.plot(x, y, ',', color="black")

			for i, (mu0, sigma0, energy, sigma_energy) in enumerate(meta["peaks"], 1):
				fit = local_fit(x,y,errors, mu=mu0, sigma=sigma0)
				fit.plot(fit.mu-5*fit.sigma, fit.mu+5*fit.sigma, 1000, zorder=10000, color="red")
				mus.append((fit.mu, fit.error_mu))

				if meta["element"] in ("Na", "Cs"):
					E_gamma = 0.09863 * fit.mu - 19.056
					m_electron = 511 # keV
					E_compton = E_gamma / (1 + m_electron / (2*E_gamma))

					c_compton = (E_compton+19.056) / 0.09863

					plt.axvline(c_compton, color="red", linestyle="--")

				peaktable.row((meta["element"], len(meta["peaks"])),
					formatQuantityLatex(fit.A, fit.error_A, parenthesis=False, suffix=False),
					formatQuantityLatex(fit.mu, fit.error_mu, parenthesis=False, suffix=False),
					formatQuantityLatex(fit.sigma, fit.error_sigma, parenthesis=False, suffix=False),
					"%.2f" % fit.chi2ndf,
					formatQuantityLatex(energy, sigma_energy, parenthesis=False, suffix=False)
				)

				widths.append((fit.sigma, fit.error_sigma))
				energies.append((energy, sigma_energy))
				print(meta["element"], i, fit)

			plt.xlabel("Kanal")
			plt.ylabel("Count")
			plt.xlim(0, 2**14)
			plt.title(meta["title"])
			plt.savefig("out/" + clean_filename(filename) + "_all." + SAVETYPE)

			N = ufloat(y.sum(), errors.sum())
			t = ufloat(tka.real_time, 1/math.sqrt(12))
			m = N/t # Hz/Bq
			#r = ufloat(91.5, 0.01) # mm
			#A = ufloat(*meta["activity"])*1000 # Bq
			#I_g = 1
			#d = ufloat(20, 10)
			#F_D = r / d * ufloat(7.45, 0.05) * ufloat(80.75, 0.05) # mm^2

			#efficiency = 4*math.pi * r**2 * m / (F_D * A * I_g)
			ratetable.row(meta["element"], formatUFloatLatex(N), formatUFloatLatex(t, unit="s"), formatUFloatLatex(m, unit="1/s"))

	mus, error_mus = np.array(mus).T
	widths, error_widths = np.array(widths).T
	energies, error_energies = np.array(energies).T

	# Kalibration
	plt.clf()
	fit = Fit(LINEAR)
	fit.slope = 0.1
	fit.offset = -20
	for i in range(5):
		errors = fit.combine_errors(mus, error_mus, error_energies)
		fit.fit(mus, energies, errors)

	plt.errorbar(mus, energies, xerr=error_mus, yerr=error_energies, fmt=',', color="black")
	fit.plot(box='tl', units={"slope": "eV / Channel", "offset": "eV"}, factors={"slope": 1000, "offset": 1000}, color="red")
	plt.xlabel("Kanal")
	plt.ylabel("Energie / keV")
	plt.savefig("out/calib_fit." + SAVETYPE)

	plt.clf()
	errs = fit.combine_errors(mus, xerr=error_mus, yerr=error_energies)
	fit.plot_residual(mus, energies, errs, box='tl', fmt=',', color="black")
	plt.xlabel("Kanal")
	plt.ylabel("Energie / keV")
	plt.savefig("out/calibresiduum." + SAVETYPE)

	s = formatQuantityLatex(fit.slope*1000, fit.error_slope*1000, unit="eV / Kanal", math=False)
	o = formatQuantityLatex(fit.offset*1000, fit.error_offset*1000, unit="eV", math=False)

	with open('out/calib.tex', 'w') as file:
		file.write(r'Die in diesem Versuch verwendeten Einstellungen f\"uhren zu einer Einteilung von \[ ' +s+ r' \] wobei der erste Kanal die Energie \[ '+o+r' \] besitzt.')

	# Energieauflösung
	error_widths = np.sqrt(np.power(fit.error_slope * widths, 2) + np.power(error_widths * fit.slope, 2))
	widths = widths * fit.slope

	error_energies = np.sqrt(np.power(fit.error_offset,2) + np.power(mus*fit.error_slope, 2) + np.power(fit.slope*error_mus,2))
	energies = fit.slope * mus  + fit.offset

	X = energies
	Y = widths
	SX = error_energies
	SY = error_widths

	func = lambda E, a, b: np.sqrt(np.power(a*E,2) + np.power(b,2)*E)
	fit = Fit(func)
	fit.a = 0.01
	fit.b = 0.01

	for _ in range(10):
		err = fit.combine_errors(X, SX, SY)
		fit.fit(X, Y, err)

	plt.clf()
	plt.errorbar(X, Y, xerr=SX, yerr=SY, fmt=',', color="black")
	fit.plot(box='br', units={'b': r'\sqrt{\textrm{keV}}'}, color="red")
	plt.xlabel(r"$ E $ / keV")
	plt.ylabel(r"$ \Delta E / keV$")
	plt.title(r'Energieauflösung: Fit zu $ \Delta E = \sqrt{\left(a E\right)^2 + b^2 E} $')
	plt.savefig("out/energyresolution_fit." + SAVETYPE)

	print("a =", formatQuantity(fit.a, fit.error_a))
	print("b =", formatQuantity(fit.b, fit.error_b))

	plt.clf()
	err = fit.combine_errors(X, SX, SY)
	fit.plot_residual(X, Y, err, box='tr', fmt=",", color="black")
	plt.xlabel(r"$ E $ / keV")
	plt.ylabel(r"$ \Delta E / keV$")
	plt.title("Energieauflösung: Residuen")
	plt.savefig("out/energyresolution_residual." + SAVETYPE)

	with LatexTable("out/energyresolution.tex") as table:
		table.header("Parameter", "Wert")
		table.row("$a$", formatQuantityLatex(fit.a, fit.error_a))
		table.row("$b$", formatQuantityLatex(fit.b, fit.error_b, unit="\sqrt{keV}"))

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

	plt.clf()
	plt.errorbar(time, count, yerr=error_count, fmt='s', color="black")

	func = lambda x, offset, mu, sigma, A: offset + GAUSS(x, mu, sigma, A)

	fit = Fit(func)
	fit.offset = 0
	fit.mu = -45
	fit.sigma = 60
	fit.fit(time, count, error_count)
	fit.plot(box="tr", units={"mu": "ns", "sigma": "ns"}, color="red")
	print(fit)

	plt.xlabel(r"$\Delta T / \textrm{ns}$")
	plt.ylabel("Count")
	plt.savefig("out/delay_fit." + SAVETYPE)

	plt.clf()
	fit.plot_residual(time, count, error_count, fmt="s", box="br", color="black")
	plt.xlabel(r"$\Delta T / \textrm{ns}$")
	plt.ylabel("Count")
	plt.savefig("out/delay_residual." + SAVETYPE)

def plot_coincidence_na():
	for spalte in (1,2):
		columns = np.loadtxt("data/natrium_winkel.txt", unpack=True)
		angle = columns[0]
		count = columns[spalte]
		error_angle = 0.3
		error_count = np.sqrt(count + 1)

		plt.clf()
		plt.errorbar(angle, count, xerr=error_angle, yerr=error_count, fmt=',', color="black")

		fit = Fit(GAUSS)
		fit.mu = 0
		fit.sigma = 10
		for i in range(5):
			errors = fit.combine_errors(angle, error_angle, error_count)
			fit.fit(angle, count, errors)
		fit.plot(fit.mu-5*fit.sigma, fit.mu+5*fit.sigma, box="bl", color="red")

		plt.ylim(0, plt.ylim()[1])
		plt.xlabel("Winkel / Grad")
		plt.ylabel("Count")
		plt.savefig("out/coincidence_na_%d_fit." % spalte + SAVETYPE)

		plt.clf()
		fit.plot_residual(angle, count, errors, fmt="s", box="br", color="black")
		plt.xlabel("Winkel / Grad")
		plt.ylabel("Count")
		plt.savefig("out/coincidence_na_%d_residual." % spalte + SAVETYPE)

def plot_coincidence_co():
	angle, count = np.loadtxt("data/winkel_cobalt3.txt", unpack=True)

	error_angle = 0.3
	error_count = np.sqrt(count + 1)

	plt.clf()
	plt.errorbar(angle, count, xerr=error_angle, yerr=error_count, fmt='s', color="black")


	W = lambda theta, A, a2, a4: A*(1 + a2*np.power(np.cos(np.radians(theta)), 2) + a4*np.power(np.cos(np.radians(theta)), 4))

	fit = Fit(W)
	fit.A = count.mean()
	fit.a2 = 0
	fit.a4 = 0
	for i in range(5):
		errors = fit.combine_errors(angle, error_angle, error_count)
		fit.fit(angle, count, errors)
	fit.plot(-90, 90, box="br", color="red")
	plt.xlim(-90, 90)
	plt.xlabel("Winkel / Grad")
	plt.ylabel("Count")
	plt.savefig("out/coincidence_co_fit." + SAVETYPE)

	plt.clf()
	fit.plot_residual(angle, count, errors, fmt="s", box="br", color="black")
	plt.xlabel("Winkel / Grad")
	plt.ylabel("Count")
	plt.savefig("out/coincidence_co_residual." + SAVETYPE)


if __name__=="__main__":
	plot_energy_spectra()
	plot_delay_spectrum()
	plot_coincidence_na()
	plot_coincidence_co()
