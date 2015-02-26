import math
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from mca import McaFile
from dft import dft, smooth
from funcs import GAUSS
from fits import linear_fit, local_gauss_fit, local_double_gauss_fit

from manuell import *

PATH = os.path.dirname(os.path.realpath(__file__))

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

		plt.errorbar(self.channels, self.data, yerr=errors_no_0, fmt='s', color="red")
		plt.xlabel("Channel")
		plt.xlim(np.min(self.channels), np.max(self.channels))
		plt.ylabel("Count")
		plt.title(self.title)

	def find_peak(self, mu0, sigma0, plot=False):
		fit = local_gauss_fit(self.data, mu0, sigma0, A=100, errors=self.errors)
		if plot:
			fit.plot(mu0-4*sigma0, mu0+4*sigma0, N=500, color="black")
		return fit

	def find_double_peak(self, mu1, sigma1, mu2, sigma2, plot=False):
		fit = local_double_gauss_fit(self.data, mu1, mu2, sigma1, sigma2, A1=100, A2=100, errors=self.errors)
		if plot:
			lower = min(mu1-4*sigma1, mu2-4*sigma2)
			upper = max(mu1+4*sigma1, mu2+4*sigma2)
			fit.plot(lower, upper, N=500, color="black")
		return fit


def kalibrierung(plot=True):
	print("##### CALIBRATION #####")
	channels = []
	channel_errors = []
	energies = []

	for filename, meta in kalib.items():
		print("="*10, filename, "="*10)
		experiment = Experiment(PATH + "/data/" + filename, title=filename)
		if plot: experiment.errorplot()
		for i, peak in enumerate(sorted(meta["peaks"], key=lambda peak: peak[0]), 1):
			if len(peak) == 3:
				mu0, sigma0, energy = peak
				fit = experiment.find_peak(mu0, sigma0, plot=plot)
				channels.append(fit.mu)
				channel_errors.append(fit.sigma_mu)
				energies.append(energy)
				print("PEAK %d: channel %.1f, width %.1f, height: %.2f, chisq: %.2f, energy: %.2f keV" % (i, fit.mu, fit.sigma, fit.A, fit.chisqndof, energy))
			elif len(peak) == 6:
				mu1, sigma1, energy1, mu2, sigma2, energy2 = peak
				fit = experiment.find_double_peak(mu1, sigma1, mu2, sigma2, plot=plot)
				channels.append(fit.mu1)
				channels.append(fit.mu2)
				channel_errors.append(fit.sigma_mu1)
				channel_errors.append(fit.sigma_mu2)
				energies.append(energy1)
				energies.append(energy2)
				print("PEAK %da: channel %.1f, width %.1f, height: %.2f, chisq: %.2f, energy: %.2f keV" % (i, fit.mu1, fit.sigma1, fit.A1, fit.chisqndof, energy1))
				print("PEAK %db: channel %.1f, width %.1f, height: %.2f, chisq: %.2f, energy: %.2f keV" % (i, fit.mu2, fit.sigma2, fit.A2, fit.chisqndof, energy2))

		if plot: plt.show()

	X = np.array(energies)
	SX = 0.01/math.sqrt(12) #keV
	Y = np.array(channels)
	SY = np.array(channel_errors)

	plt.errorbar(X, Y, xerr=SX, yerr=SY, fmt='s')

	fit = linear_fit(X, Y, xerr=SX, yerr=SY)
	if plot:
		fit.plot(np.min(X), np.max(X))
		plt.xlabel("Energy / keV")
		plt.ylabel("Channel")
		plt.show()

		err = fit.combine_errors(X, SX, SY)
		fit.plot_residuums(X, Y, err, fmt="s")
		plt.xlabel("Energy / keV")
		plt.ylabel("Channel")
		plt.show()

	slope, offset = 1/fit.slope, -fit.offset/fit.slope

	print("RESULT: linear fit with chisq: %.2f" % (fit.chisqndof))
	print("RESULT: energy per bin: %.1f eV" % (slope*1000))
	print("RESULT: energy of 0th bin: %.1f eV" % (offset*1000))
	print()
	return slope, offset


def auswertung(plot=True):
	slope, offset = kalibrierung(plot=plot)

	print("##### UNKNOWN #####")
	for filename, meta in data.items():
		print("="*10, filename, "="*10)
		experiment = Experiment(PATH + "/data/" + filename, title=filename, energy_slope=slope, energy_offset=offset)

		experiment.subtract_empty(PATH + "/data/G20_Leer.mca", 0.5)
		if plot: experiment.errorplot()

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
			plt.show()
	print()

if __name__=="__main__":
	# for filename in kalib:
	# 	experiment = Experiment("data/" + filename)
	# 	experiment.errorplot()
	# 	plt.show()
	auswertung(True)
	#kalibrierung()
	#testsignal()
