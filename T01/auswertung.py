import os
import math

from pylab import rcParams
rcParams['savefig.dpi'] = 200
#rcParams['text.usetex'] = True
rcParams['text.latex.unicode'] = True
#rcParams['text.latex.preamble']=[r'\usepackage{lmodern}']
#rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'Computer Modern'
rcParams['font.size'] = 18
rcParams['savefig.bbox'] = 'tight'
rcParams['axes.grid'] = True

import numpy as np
import matplotlib.pyplot as plt

from uncertainties import wrap as uwrap
from uncertainties import ufloat, umath

from funcs import *
from fits import Fit, info_box
from number import *
from tka import *

from scipy.interpolate import interp1d

os.chdir(os.path.dirname(os.path.realpath(__file__)))

SAVETYPE = "png"

def plot_theo_alpha():
	peaks = np.array([4.871, 5.307, 5.590, 6.115, 7.69])
	widths = 0.07 * np.sqrt(peaks)

	x = np.linspace(0,10, 1000)
	y = 0
	y2 = 0
	for peak, width in zip(peaks, widths):
		y += GAUSS(x, peak, width)
		y2 += GAUSS(x, peak, 0.1*width)

	plt.clf()
	plt.minorticks_on()
	plt.plot(x,y)
	plt.plot(x,y2)
	plt.xlabel("Energie / MeV")
	plt.savefig("out/theo." + SAVETYPE)

def plot_radium_closest():
	tka = TkaFile("data/RadiumNaechsterAbstand.tka")
	channels = np.arange(1,len(tka)+1)

	plt.clf()
	plt.minorticks_on()
	plt.plot(channels, tka.data, ".", alpha=0.7)
	plt.xlabel("Kanal")
	plt.ylabel("Anzahl")
	plt.savefig("out/radium_close." + SAVETYPE)
	plt.xlim(0,2000)
	plt.ylim(0,2000)
	plt.savefig("out/radium_close2." + SAVETYPE)
	plt.show()

def plot_radium_empty():
	tka = TkaFile("data/RadiumLeermessung.tka")
	channels = np.arange(1,len(tka)+1)

	plt.clf()
	plt.minorticks_on()
	plt.plot(channels, tka.data, ".", alpha=0.7)
	plt.xlabel("Kanal")
	plt.ylabel("Anzahl")
	plt.savefig("out/radium_empty." + SAVETYPE)


def plot_radium_peak_distances():
	abstand = {1: 35.3, 2: 36.2, 3: 36.9, 4: 39}
	for i in range(1,5):
		plt.clf()
		plt.minorticks_on()
		for pos in ("vor", "hinter", "drauf"):
			filename = "data/RadiumPeak%d%s.TKA" % (i, pos)
			tka = TkaFile(filename)
			channels = np.arange(1, len(tka)+1)
			symbol, legend = {
				"vor": (".", "1 mm näher"),
				"drauf": (".", "exakt %.1f cm" % abstand[i]),
				"hinter": (".", "1 mm weiter"),
			}[pos]
			plt.plot(channels, tka.data, symbol, label=legend, alpha=0.7)
		plt.xlabel("Kanal")
		plt.ylabel("Anzahl")
		plt.xlim(0,2000)
		plt.ylim(0,{
			1: 2000,
			2: 1000,
			3: 500,
			4: 250,
		}[i])
		plt.legend()
		plt.savefig("out/radium_peak%d." % (i, ) + SAVETYPE)

def plot_radium_distance():
	for distance in ("34_9", "35_0", "35_5", "36_0", "36_5", "37_0", "37_5", "38_0", "38_5", "39_0", "39_5"):
		filename = "data/Radium%s.TKA" % distance
		tka = TkaFile(filename)
		channels = np.arange(1, len(tka)+1)
		plt.clf()
		plt.minorticks_on()
		plt.plot(channels, tka.data, ".", alpha=0.7)
		plt.xlabel("Kanal")
		plt.ylabel("Anzahl")
		plt.xlim(0,2000)
		plt.ylim(0,2000)
		info_box("Distance: " + distance.replace("_", ".") + "cm", location="tr")
		plt.savefig("out/radium_distance_%s." % distance + SAVETYPE)

def fit_radium_calib():
	energies = [3.953, 5.332, 7.02] # MeV
	error_energies = np.array([0.001, 0.001, 0.01]) / math.sqrt(12)
	channels = [240, 815, 1495]
	error_channels = [10, 10, 8]
	fit = Fit(LINEAR)
	fit.set_data(xdata=channels, xerrors=error_channels, ydata=energies, yerrors=error_energies)
	fit.set_labels(xlabel="Kanal", ylabel="Energie / MeV")
	fit.iterative_fit(5)
	plt.clf()
	plt.minorticks_on()
	fit.plot(box="tl", units={"slope": "keV/Kanal", "offset": "MeV"}, factors={"slope": 1000})
	plt.savefig("out/radium_calib_fit." + SAVETYPE)
	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="bl")
	plt.savefig("out/radium_calib_residual." + SAVETYPE)

def radium_calib_2():
	E, proj_range = np.loadtxt("ranges.txt", unpack=True)
	real_range = proj_range / 1.184e-3

	energy_from_range = interp1d(real_range, E, kind='cubic')

	energies = []
	error_energies = []
	for distance in (0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5):
		distance = ufloat(35 + distance, 0.05)
		total_distance = (distance-35) + ufloat(3.17, 0.17)
		for peak in (3.52, 4.99, 7.18):
			peak = ufloat(peak, 0.01)

			rest = peak - total_distance
			try:
				upper = energy_from_range(rest.n + rest.s)
				lower = energy_from_range(rest.n - rest.s)
				value = energy_from_range(rest.n)
			except ValueError as e:
				#print(e)
				#print("orig distance:", peak, "cm  ", "rest distance:", rest, "cm  ")
				pass
			else:
				error = max(upper-value, value-lower)
				energy = ufloat(value, error)
				#print("orig distance:", peak, "cm  ", "total distance:", total_distance, "cm  ", "rest distance:", rest, "cm  ", "rest energy:", energy, "MeV")
				print("$({:L}) \\unit{{cm}}$ & $({:L}) \\unit{{cm}}$ & $({:L}) \\unit{{cm}}$ & $({:L}) \\unit{{cm}}$ & $({:L}) \\unit{{MeV}}$ \\\\".format(peak, distance, total_distance, rest, energy))

				energies.append(value)
				error_energies.append(error)


	# PART2 !!!
	channels = [220, 780, 1460, 620, 1340, 450, 1260, 150, 1170, 990, 820, 620, 390]
	error_channels = [20, 20, 20, 20, 20, 20, 20, 100, 20, 20, 20, 20, 30]
	energies = energies[:len(channels)]
	error_energies = error_energies[:len(channels)]

	for energy, error_energy, channel, error_channel in zip(energies, error_energies, channels, error_channels):
		e = ufloat(energy, error_energy)
		c = ufloat(channel, error_channel)
		print("$({:L}) \\unit{{MeV}}$ & ${:L}$ \\\\".format(e,c))

	fit = Fit(LINEAR)
	fit.set_data(xdata=channels, xerrors=error_channels, ydata=energies, yerrors=error_energies)
	fit.set_labels(xlabel="Kanal", ylabel="Energie / MeV")
	fit.iterative_fit(5)
	plt.clf()
	plt.minorticks_on()
	fit.plot(box="tl", units={"slope": "keV/Kanal", "offset": "MeV"}, factors={"slope": 1000})
	plt.savefig("out/radium_calib2_fit." + SAVETYPE)
	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tr")
	plt.savefig("out/radium_calib2_residual." + SAVETYPE)

def calc_radium_range():
	peak = {
		35:	(1100, 1800),
		35.5: (1000, 1700),
		36: (900, 1600),
		36.5: (800, 1500),
		37: (600, 1400),
		37.5: (500, 1200),
		38: (300, 1000),
		38.5: (100, 700)
	}

	distances = []
	error_distances = []
	integrals = []
	error_integrals = []

	for distanceStr in ("35_0", "35_5", "36_0", "36_5", "37_0", "37_5", "38_0", "38_5"):
		filename = "data/Radium%s.TKA" % distanceStr

		distance = ufloat(float(distanceStr.replace("_", ".")), 0.05)

		tka = TkaFile(filename)
		lower, upper = peak[distance.n]

		data = tka.data[lower:upper]

		distance = distance - 35 + ufloat(3.17, 0.17)

		integral = ufloat(data.sum(), math.sqrt((data+1).sum()))
		integral_fixed = integral * distance**2

		print("$({:L}) \\unit{{cm}}$ & {:d} & {:d} &  ${:dL}$ &  ${:dL} \\unit{{cm^2}}$ \\\\".format(distance, lower, upper, integral, integral_fixed))

		distances.append(distance.n)
		error_distances.append(distance.s)
		integrals.append(integral_fixed.n)
		error_integrals.append(integral_fixed.s)

	distances = np.array(distances)
	integrals = np.array(integrals)
	error_distances = np.array(error_distances)
	error_integrals = np.array(error_integrals)

	plt.clf()
	plt.minorticks_on()
	plt.errorbar(distances, integrals, xerr=error_distances, yerr=error_integrals, fmt="s")

	#plt.plot(distances, integrals, 's')
	plt.xlabel("Distanz / cm")
	plt.ylabel("korrigierte Summe / cm^2")
	plt.ylim(0, plt.ylim()[1])
	plt.savefig("out/radium_range." + SAVETYPE)
	#plt.show()

	plt.clf()
	plt.minorticks_on()
	scale = 1/integrals[0:3].mean()

	fit = Fit(LINEAR)
	fit.set_params(offset=1, slope=-0.5)
	fit.set_data(xdata=distances[1:], ydata=integrals[1:]*scale, xerrors=error_distances[1:], yerrors=error_integrals[1:]*scale)
	fit.set_labels(xlabel="Distanz / cm", ylabel="Intensität")
	fit.iterative_fit(5)
	fit.plot(plot_data=True, plot_fit=True, box="bl", units={"slope":"1/cm"})
	plt.ylim(0, plt.ylim()[1])
	plt.savefig("out/radium_range_fit." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tr")
	plt.savefig("out/radium_range_residual." + SAVETYPE)

	r = (0.5-fit.uvalue("offset"))/fit.uvalue("slope")
	print("Range:", r, "cm")


def plot_ionisation_raw():
	distances, voltages, error_voltages = np.loadtxt("data/radium_ionisationskammer.txt", unpack=True)

	# discard 38.80cm
	distances = distances[:-1]
	voltages = voltages[:-1]
	error_voltages = error_voltages[:-1]

	for distance, voltage, error_voltage in zip(distances, voltages, error_voltages):
		distance = ufloat(distance, 0.05)
		voltage = ufloat(voltage, error_voltage)

		print("$({:L}) \\unit{{cm}}$ & $({:L}) \\unit{{mV}}$ \\\\".format(distance, voltage))

	plt.clf()
	plt.minorticks_on()
	plt.errorbar(distances, voltages, xerr=0.05, yerr=error_voltages, fmt=',')
	plt.xlabel("Distanz / cm")
	plt.ylabel("Spannung / mV")
	plt.savefig("out/radium_ionisation_raw." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()

	currents = []
	error_currents = []
	for voltage, error_voltage in zip(voltages, error_voltages):
		current = ufloat(voltage, error_voltage) #+ ufloat(5, 3)
		currents.append(current.n)
		error_currents.append(current.s)

	currents = np.array(currents)
	error_currents = np.array(error_currents)

	distances -= 39

	plt.clf()
	plt.minorticks_on()
	plt.errorbar(distances, currents, xerr=0.05, yerr=error_currents, fmt=',')
	plt.xlabel("Distanz / cm")
	plt.ylabel("Strom / nA")
	plt.savefig("out/radium_ionisation." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()

	x = distances[1:] + np.abs(np.diff(distances)/2)
	y = - np.diff(currents) / np.diff(distances)
	plt.plot(x, y, 's-')

	plt.xlabel("Distanz / cm")
	plt.ylabel("- Änderung des Stromes / nA / cm")
	plt.savefig("out/radium_ionisation_diff." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	fit = Fit(LINEAR)
	fit.set_data(xdata=(7.42, 6.92, 6.42), ydata=currents[9:12], xerrors=0.25, yerrors=error_currents[9:12])
	fit.set_labels(xlabel="Distanz / cm", ylabel="Strom / nA")
	fit.iterative_fit(5)
	fit.plot(box="tr", units={"slope": "nA/cm", "offset": "nA"})
	plt.savefig("out/radium_ionisation_po_fit." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tl")
	plt.savefig("out/radium_ionisation_po_residual." + SAVETYPE)

	middle = (ufloat(currents[9], error_currents[9])+ufloat(currents[11], error_currents[11])) / 2
	print("Middle:",middle, "nA")
	r = (middle - fit.uvalue("offset"))/fit.uvalue("slope")
	print("Range:",r,"cm")

def plot_caesium_spectra():
	for scintillator in ("Organisch", "Anorganisch"):
		print("===", scintillator, "===")
		tka = TkaFile("data/Caesium_10min_" + scintillator + ".TKA")
		channels = np.arange(1,len(tka)+1)
		plt.clf()
		plt.minorticks_on()
		plt.plot(channels, tka.data, '.')
		plt.xlabel("Kanal")
		plt.ylabel("Anzahl")
		plt.savefig("out/caesium_" + scintillator + "_raw." + SAVETYPE)
		plt.ylim(0,2500)
		plt.savefig("out/caesium_" + scintillator + "." + SAVETYPE)
		plt.show()

		print("Summe:", tka.data.sum())

		#
		# fit = Fit(GAUSS)
		#
		# if scintillator == "Organisch":
		# 	lower = 100
		# 	upper = 5000
		# elif scintillator == "Anorganisch":
		# 	lower = 100
		# 	upper = 500
		#
		# print("Obere Summe:", tka.data[upper:].sum())
		#
		# fit.set_data(xdata=channels[lower:upper], ydata=tka.data[lower:upper], yerrors=np.sqrt(tka.data[lower:upper]+1))
		# fit.set_labels(xlabel="Kanal", ylabel="Anzahl")
		# fit.iterative_fit(5)
		#
		# data = tka.data - fit.eval(channels)
		# plt.clf()
		# plt.minorticks_on()
		# plt.plot(channels, data, '.')
		# plt.xlabel("Kanal")
		# plt.ylabel("Anzahl")
		# plt.ylim(0,2500)
		# plt.savefig("out/caesium_" + scintillator + "_corrected." + SAVETYPE)

def plot_caesium_absorber():
	depth, events, slabs = np.loadtxt("data/caesium_bleiabschirmung.txt", unpack=True)
	error_depth = np.sqrt(slabs) * 0.1
	error_events = np.sqrt(events + 1)

	plt.clf()
	plt.minorticks_on()
	plt.errorbar(depth, events, xerr=error_depth, yerr=error_events, fmt=',')
	plt.xlabel("Tiefe / mm")
	plt.ylabel("ln(Anzahl)")
	plt.savefig("out/caesium_absorber_nolog." + SAVETYPE)

	ln_events = np.log(events)
	error_ln_events = np.abs(1/events)*error_events

	depth /= 10 # mm=>cm
	error_depth /= 10

	plt.clf()
	plt.minorticks_on()

	func = lambda x, mu, lnI0: -x*mu + lnI0

	fit = Fit(func)
	fit.set_data(xdata=depth, ydata=ln_events, xerrors=error_depth, yerrors=error_ln_events)
	fit.set_labels(xlabel="Tiefe / cm", ylabel="ln(Anzahl)")
	fit.iterative_fit(5)
	fit.plot(box="tr", units={"mu": "1/cm"})
	plt.savefig("out/caesium_absorber_fit." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tr")
	plt.savefig("out/caesium_absorber_residual." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	fit = fit.filtered(depth <= 2)
	fit.iterative_fit(5)
	fit.plot(box="tr", units={"mu": "1/cm"})
	plt.savefig("out/caesium_absorber_fit2." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tr")
	plt.savefig("out/caesium_absorber_residual2." + SAVETYPE)

	mu = fit.uvalue("mu")
	print("Absorptionskoeffizient:", mu, "1/cm")
	mu1 = mu / ufloat(11.342, 0.001)
	print("Massenabsorptionskoeffizient:", mu1, "cm^2/g")


def caesium_absorb_nist():
	energy, absorb, _ = np.loadtxt("lead.txt", unpack=True)
	# plt.clf()
	# plt.minorticks_on()
	# plt.plot(energy, absorb)
	# plt.xlabel("Energy / MeV")
	# plt.ylabel("Absorption / cm^2/g")
	# plt.xscale("log")
	# plt.yscale("log")
	# plt.show()

	part = np.logical_and(energy > 0.1, energy < 10)
	func = interp1d(energy[part], absorb[part], kind='cubic')

	E = 0.6617
	mu1 = func(E)
	print("theor. Massenabsorptionskoeffizient:", mu1, "cm^2/g")
	mu = mu1 * 11.342
	print("theor. Absorptionskoeffizient:", mu, "1/cm")

def plot_strontium():
	m = {
		0:	0,
		5:	0.03,
		6:	0.04,
		7:	0.05,
		8:	0.08,
		9:	0.10,
		10:	0.17,
		11: 0.25,
		12:	0.37,
		13:	0.50,
		14:	0.64,
		15:	0.83,
		16:	1.01,
		17:	1.23,
		18:	1.39,
		19:	1.57,
		20:	1.88,
		21:	2.29,
		22:	2.79,
		23:	3.46,
	}

	nr, count = np.loadtxt("data/strontium_alu.txt", unpack=True)
	absorb = np.zeros_like(nr)
	for i, x in enumerate(nr):
		absorb[i] = m[x] / 10

	scale = 1/count[0]

	error_count = np.sqrt(count+ 1)
	error_absorb = 0.01 / 10


	count *= scale
	error_count *= scale

	func = lambda x, A, mu, offset: A*np.exp(-x*mu) + offset

	fit = Fit(func)
	fit.set_data(xdata=absorb, ydata=count, xerrors=error_absorb, yerrors=error_count)
	fit.set_params(A=1, mu=1, offset=0)
	fit.set_labels(xlabel="Dicke / cm", ylabel="Anteil")
	fit.iterative_fit(5)

	plt.clf()
	plt.minorticks_on()
	plt.xlim(0,.35)
	fit.plot(box="tr", units={"mu": "1/cm"})
	plt.savefig("out/strontium_fit." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tr")
	plt.savefig("out/strontium_residual." + SAVETYPE)

	rho = ufloat(2.7, 0.1)
	mu = fit.uvalue("mu")
	print("Absorptionskoeffizient:", mu, "1/cm")
	mu1 = mu / rho
	print("Massenabsorptionskoeffizient:", mu1, "cm^2/g")

	E = umath.pow(17 / mu1, 1/1.14)
	print("Maximalenergie:", E, "MeV")

	R1 = 0.412*umath.pow(E,1.265-0.0954*umath.log(E))
	R = R1 / rho

	print("Reichweite:", R1, "g/cm^2")
	print("Reichweite:", R, "cm")

if __name__=="__main__":
	plot_theo_alpha()
	plot_radium_closest()
	plot_radium_peak_distances()
	plot_radium_distance()
	fit_radium_calib()
	radium_calib_2()
	calc_radium_range()
	plot_ionisation_raw()
	plot_caesium_spectra()
	plot_caesium_absorber()
	caesium_absorb_nist()
	plot_strontium()
