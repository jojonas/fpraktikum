import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from mca import McaFile
from dft import dft, smooth
from funcs import GAUSS
from peaks import gauss_fit, chisqNdof
from linreg import linreg2

from manuell import *

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
	LEER = McaFile("data/G20_Leer.mca")
	LEER_FACTOR = 0.5

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
		mu, sigma, A = gauss_fit(self.data, mu0, sigma0, 100, errors=self.errors)
		if plot:
			self.plot_peak(mu, sigma, A)
		return mu, sigma, A

	def chisq_peak(self, mu, sigma, A):
		lower = max(int(mu-sigma), 0)
		upper = min(int(mu+sigma), self.count)
		xdata = self.channels[lower:upper]
		ydata = self.data[lower:upper]
		yerrors = self.errors[lower:upper]
		return chisqNdof(xdata, ydata, lambda x: GAUSS(x,A,mu,sigma), yerrors, nparams=3)

	def plot_peak(self, mu, sigma, A, color="black", npoints=100):
		nx = np.linspace(mu-5*sigma, mu+5*sigma, npoints)
		ny = GAUSS(nx, A, mu, sigma)
		plt.plot(nx, ny, '-', color=color)

		arrow_y = A*np.exp(-1/2)
		opts = {
			"length_includes_head": True,
			"head_width": 4,
			"head_length": 4,
			"fc": color,
			"ec": color,
			"zorder": 100000
		}

		plt.axes().arrow(mu-sigma, arrow_y, 2*sigma, 0, **opts)
		plt.axes().arrow(mu+sigma, arrow_y, -2*sigma, 0, **opts)
		plt.axes().arrow(mu, 0, 0, A, **opts)

def kalibrierung():
	channels = []
	channel_errors = []
	energies = []

	for filename, meta in kalib.items():
		experiment = Experiment("data/" + filename, title=filename)
		experiment.errorplot()
		for mu0, sigma0, energy in meta["peaks"]:
			mu, sigma, A = experiment.find_peak(mu0, sigma0)
			experiment.plot_peak(mu, sigma, A)
			#chisq = experiment.chisq_peak(mu, sigma, A)
			channels.append(mu)
			channel_errors.append(sigma)
			energies.append(energy)
		plt.show()

	X = energies
	SX = 0.01/math.sqrt(12) #keV
	Y = channels
	SY = channel_errors

	plt.errorbar(X, Y, xerr=SX, yerr=SY, fmt='s')
	popt, perr = linreg2(X, Y, SX, SY)
	slope, offset = popt
	x = np.linspace(np.min(X), np.max(X), 100)
	plt.plot(x, slope*x+offset, '-')
	plt.xlabel("Energy / keV")
	plt.ylabel("Channel")
	plt.show()

	return 1/slope, -offset/slope


def auswertung():
	slope, offset = kalibrierung()
	#slope = 0.0133341236267
	#offset = 0.0409739312723

	print("slope =", slope)
	print("offset =", offset)

	for filename, meta in data.items():
		print("="*10, filename, "="*10)
		experiment = Experiment("data/" + filename, title=filename, energy_slope=slope, energy_offset=offset)

		experiment.subtract_empty("data/G20_Leer.mca", 0.5)
		experiment.errorplot()

		for mu0, sigma0 in meta["peaks"]:
			mu, sigma, A = experiment.find_peak(mu0, sigma0)
			experiment.plot_peak(mu, sigma, A)
			chisq = experiment.chisq_peak(mu, sigma, A)

			mu = experiment.channel2energy(mu)
			sigma = experiment.channel2energy(sigma)
			print("Peak at %.3f keV with std %.3f keV and height %d... Chi^2/ndof: %.2f" % (mu, sigma, A, chisq))

		experiment.set_energy_labels()

		plt.show()

if __name__=="__main__":
	# for filename in kalib:
	# 	experiment = Experiment("data/" + filename)
	# 	experiment.errorplot()
	# 	plt.show()
	auswertung()
	#kalibrierung()
	#testsignal()
