import math
import os

from pylab import rcParams
rcParams['savefig.dpi'] = 100

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from mca import McaFile
from dft import dft, smooth
from funcs import *
from fits import Fit, linear_fit, local_gauss_fit

from manuell import *

PATH = os.path.dirname(os.path.realpath(__file__))

SAVETYPE = "png"

def testsignal():
	N = 4096
	signal = np.zeros(N)
	x = np.arange(N)
	signal += GAUSS(x, 700, 600, 20)
	signal += GAUSS(x, 150, 500, 10)
	signal += GAUSS(x, 150, 660, 10)
	f, fy = dft(signal)
	plt.plot(f, np.abs(fy), 's')
	plt.show()

class Experiment:
	def __init__(self, filename=None, title=None, energy_slope=1, energy_offset=0):
		if filename:
			file = McaFile(filename)
			self._data = file.data
		else:
			self._data = None

		if title:
			self.title = title
		else:
			self.title = filename

		self._errors = np.sqrt(self._data)
		self._energy_slope = energy_slope
		self._energy_offset = energy_offset

	@property
	def count(self):
		return len(self.data)

	@property
	def channels(self):
		return np.arange(self.count)

	@property
	def data(self):
		return self._data

	@property
	def errors(self):
		return self._errors

	@property
	def energies(self):
		return self.channel2energy(self.channels)

	def channel2energy(self, channel):
		return self._energy_slope*channel + self._energy_offset

	def channelwidth2energywidth(self, width):
		return self._energy_slope*width

	def set_energy_labels(self):
		def func(x, pos):
			return "%.2f" % self.channel2energy(x)
		formatter = FuncFormatter(func)
		plt.gca().xaxis.set_major_formatter(formatter)
		plt.xlabel("Energy / keV")

	def subtract_empty(self, filename, factor=1.0):
		file = McaFile(filename)
		self._data -= factor*file.data
		self._data[self._data<0] = 0
		self._errors = np.sqrt(np.power(self._errors, 2) + factor**2 * file.data)

	def smooth(self, fmax):
		self._data = smooth(self._data, fmax)

	def errorplot(self):
		errors_lower = self.data - np.maximum(self.data-self.errors, 0)
		errors_no_0 = np.vstack((errors_lower, self.errors))

		plt.errorbar(self.channels, self.data, yerr=errors_no_0, fmt=',', color="red")
		# plt.plot(self.channels, self.data, color="red")
		# plt.plot(self.channels, self.data + self.errors, color="red", alpha=0.5)
		# plt.plot(self.channels, self.data - errors_lower, color="red", alpha=0.5)
		plt.xlabel("Channel")
		plt.xlim(np.min(self.channels), np.max(self.channels))
		plt.ylabel("Count")
		plt.title(self.title)

	def find_peak(self, mu0, sigma0, plot=False):
		fit = local_gauss_fit(self.data, mu0, sigma0, A=100, errors=self.errors)
		if plot:
			fit.plot(fit.mu-4*fit.sigma, fit.mu+4*fit.sigma, N=500, color="black")
		return fit

def kalibration(plot=True):
	print("##### CALIBRATION #####")
	channels = []
	channel_errors = []
	energies = []

	for filename, meta in kalib.items():
		print("="*10, filename, "="*10)
		experiment = Experiment(PATH + "/data/" + filename, title=filename)
		if plot:
			plt.clf()
			experiment.errorplot()
		for i, (mu0, sigma0, energy) in enumerate(sorted(meta["peaks"], key=lambda peak: peak[0]), 1):
			fit = experiment.find_peak(mu0, sigma0, plot=plot)
			channels.append(fit.mu)
			channel_errors.append(fit.sigma_mu)
			energies.append(energy)
			print("PEAK %d: channel %.1f, width %.1f, height: %.2f, chisq: %.2f, energy: %.2f keV" % (i, fit.mu, fit.sigma, fit.A, fit.chisqndof, energy))
		if plot:
			plt.savefig("out/" + filename + "_all." + SAVETYPE)
			#plt.show()

	X = np.array(energies)
	SX = 0.1/math.sqrt(12) #keV
	Y = np.array(channels)
	SY = np.array(channel_errors)

	fit = linear_fit(X, Y, xerr=SX, yerr=SY)
	if plot:
		plt.clf()
		plt.errorbar(X, Y, xerr=SX, yerr=SY, fmt=',')
		fit.plot(np.min(X), np.max(X))
		plt.xlabel("Energy / keV")
		plt.ylabel("Channel")
		plt.savefig("out/calib_fit." + SAVETYPE)
		#plt.show()

		plt.clf()
		err = fit.combine_errors(X, SX, SY)
		fit.plot_residuums(X, Y, err, fmt=",")
		plt.xlabel("Energy / keV")
		plt.ylabel("Channel")
		plt.savefig("out/calib_residuum." + SAVETYPE)
		#plt.show()

	slope, offset = 1/fit.slope, -fit.offset/fit.slope

	print("RESULT: linear fit with chisq: %.2f" % (fit.chisqndof))
	print("RESULT: energy per bin: %.1f eV" % (slope*1000))
	print("RESULT: energy of 0th bin: %.1f eV" % (offset*1000))
	print()
	return slope, offset


def auswertung(slope, offset, plot=True):
	print("##### UNKNOWN #####")
	for filename, meta in data.items():
		print("="*10, filename, "="*10)
		experiment = Experiment(PATH + "/data/" + filename, title=filename, energy_slope=slope, energy_offset=offset)
		experiment.subtract_empty(PATH + "/data/G20_Leer.mca", 0.5)
		#experiment.smooth(0.1)
		if plot:
			plt.clf()
			experiment.errorplot()

		for i, (mu0, sigma0) in enumerate(sorted(meta["peaks"], key=lambda peak: peak[0]), 1):
			try:
				fit = experiment.find_peak(mu0, sigma0, plot=plot)
			except RuntimeError:
				print("PEAK %d: NOT FOUND!" % i)
			else:
				energy = experiment.channel2energy(fit.mu)
				width = experiment.channel2energy(fit.sigma)
				print("PEAK %d: %.3f keV, width %.3f keV, height %d, chisq %.2f" % (i, energy, width, fit.A, fit.chisqndof))

		if plot:
			experiment.set_energy_labels()
			plt.savefig("out/" + filename + "_all." + SAVETYPE)
			#plt.show()
	print()

def energieaufloesung(slope, offset, plot=True):
	energies = []
	widths = []
	sigma_energies = []
	sigma_widths = []

	for filename, meta in kalib.items():
		experiment = Experiment(PATH + "/data/" + filename, title=filename, energy_slope=slope, energy_offset=offset)
		for peak in meta["peaks"]:
			mu0, sigma0 = peak[0], peak[1]
			fit = experiment.find_peak(mu0, sigma0, plot=False)
			energies.append(experiment.channel2energy(fit.mu))
			widths.append(experiment.channelwidth2energywidth(fit.sigma))
			sigma_energies.append(experiment.channelwidth2energywidth(fit.sigma_mu))
			sigma_widths.append(experiment.channelwidth2energywidth(fit.sigma_sigma))

	for filename, meta in data.items():
		experiment = Experiment(PATH + "/data/" + filename, title=filename, energy_slope=slope, energy_offset=offset)
		experiment.subtract_empty(PATH + "/data/G20_Leer.mca", 0.5)
		for peak in meta["peaks"]:
			mu0, sigma0 = peak[0], peak[1]
			fit = experiment.find_peak(mu0, sigma0, plot=False)
			energies.append(experiment.channel2energy(fit.mu))
			widths.append(experiment.channelwidth2energywidth(fit.sigma))
			sigma_energies.append(experiment.channelwidth2energywidth(fit.sigma_mu))
			sigma_widths.append(experiment.channelwidth2energywidth(fit.sigma_sigma))


	energies = np.array(energies)
	widths = np.array(widths)
	sigma_energies = np.array(sigma_energies)
	sigma_widths = np.array(sigma_widths)

	# X = energies
	# Y = np.power(widths/energies, -2)
	#
	# SX = sigma_energies
	# SY = 2*Y * np.sqrt(np.power(sigma_widths/widths, 2) + np.power(sigma_energies/energies, 2))

	X  = energies
	Y = widths/energies
	SX = sigma_energies
	SY = Y*np.sqrt(np.power(sigma_widths/widths, 2) + np.power(sigma_energies/energies, 2))

	func = lambda x, a, b, c: a/np.sqrt(x) + b + c/x

	fit = Fit(func)
	fit.a = 1
	fit.b = 0
	for _ in range(10):
		err = fit.combine_errors(X, SX, SY)
		fit.fit(X, Y, err)

	print("RESULT: error fit with chisq: %.2f" % (fit.chisqndof))

	if plot:
		plt.clf()
		plt.errorbar(X, Y, xerr=SX, yerr=SY, fmt=',')
		fit.plot(np.min(X), np.max(X))
		plt.xlabel("Energy / keV")
		plt.ylabel("Peak width / keV^2")
		plt.savefig("out/energyresolution_fit." + SAVETYPE)

		plt.clf()
		err = fit.combine_errors(X, SX, SY)
		fit.plot_residuums(X, Y, err, fmt=",")
		plt.xlabel("Energy / keV")
		plt.ylabel("Peak width / keV^2")
		plt.savefig("out/energyresolution_residuum." + SAVETYPE)
		#plt.show()

def plot_leermessung(slope, offset):
	experiment = Experiment(PATH + "/data/G20_Leer.mca", title="Leermessung", energy_slope=slope, energy_offset=offset)
	plt.clf()
	experiment.errorplot()
	experiment.set_energy_labels()
	plt.savefig("out/leer." + SAVETYPE)


if __name__=="__main__":
	#mng = plt.get_current_fig_manager()
	#mng.window.state('zoomed')

	slope, offset = kalibration(True)
	plot_leermessung(slope, offset)
	auswertung(slope, offset, True)
	energieaufloesung(slope, offset, True)
	#testsignal()
