import math
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
from matplotlib.ticker import FuncFormatter, FixedLocator

from mca import McaFile
from dft import dft, smooth
from funcs import *
from fits import info_box, format_error, Fit, linear_fit, local_gauss_fit
from latextable import LatexTable

from number import formatNumber as _n, formatNumberPair as _p

from manuell import *

os.chdir(os.path.dirname(os.path.realpath(__file__)))

SAVETYPE = "png"

def testsignal():
	N = 4096
	signal = np.zeros(N)
	x = np.arange(N)
	signal += GAUSS(x, 601, 5.7, 2765)
	signal += GAUSS(x, 665, 6.0, 528)

	plt.clf()
	plt.plot(x, signal, '.')
	plt.xlabel("Kanal")
	plt.ylabel("simulierte Höhe")
	plt.savefig("out/testsignal." + SAVETYPE)

	f, fy = dft(signal)
	plt.clf()
	plt.plot(f, np.abs(fy), '.')
	plt.xlim(0,0.15)
	plt.axvline(0.1)
	plt.xlabel("Frequenz")
	plt.ylabel("Betrag des Fouriertransformierten")
	plt.savefig("out/testsignal_dft." + SAVETYPE)

def clean_filename(filename):
	return filename.replace(".", "_")

class Experiment:
	def __init__(self, filename=None, title=None, calibration=None):
		if filename:
			file = McaFile(filename)
			self._data = file.data
			self._mcameta = file.meta
		else:
			self._data = None

		if title:
			self.title = title
		else:
			self.title = filename

		self._errors = np.sqrt(self._data+1)
		self._calibration = calibration

	@property
	def mcameta(self):
		return self._mcameta

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
		assert self._calibration, "No calibration set."
		return self._calibration.slope*channel + self._calibration.offset

	def channelwidth2energywidth(self, width):
		assert self._calibration, "No calibration set."
		return self._calibration.slope*width

	def sigma_channel2energy(self, channel):
		return np.sqrt(np.power(channel*self._calibration.sigma_slope, 2) + np.power(self._calibration.sigma_offset, 2))

	def sigma_channelwidth2energywidth(self, width):
		return width*self._calibration.sigma_slope

	def set_energy_labels(self, stepsize=5):
		def func(x, pos):
			return "$" + _n(self.channel2energy(x), precision=3) + "$"
		formatter = FuncFormatter(func)
		plt.gca().xaxis.set_major_formatter(formatter)
		pos = (np.linspace(0, 60, 60/stepsize+1) - self._calibration.offset)/self._calibration.slope
		locator = FixedLocator(pos)
		plt.gca().xaxis.set_major_locator(locator)
		plt.xlabel("$ E $ / keV")

	def subtract_empty(self, filename, factor=1.0):
		file = McaFile(filename)
		self._data -= factor*file.data
		self._data[self._data<0] = 0
		self._errors = np.sqrt(np.power(self._errors, 2) + factor**2 * (file.data+1))

	def smooth(self, fmax):
		self._data = smooth(self._data, fmax)

	def errorplot(self):
		errors_lower = self.data - np.maximum(self.data-self.errors, 0)
		errors_no_0 = np.vstack((errors_lower, self.errors))

		plt.plot(self.channels, self.data, '.-', color="black")
		plt.fill_between(self.channels, self.data-errors_lower, self.data+self.errors, color="black", alpha=0.2)
		plt.xlabel("Channel")
		plt.xlim(np.min(self.channels), np.max(self.channels))
		plt.ylabel("Count")
		plt.title(self.title)

	def plot(self):
		plt.plot(self.channels, self.data, '-', color="black")
		plt.xlabel("Channel")
		plt.xlim(np.min(self.channels), np.max(self.channels))
		plt.ylabel("Count")
		plt.title(self.title)

	def find_peak(self, mu0, sigma0, plot=False):
		xerrs = np.ones_like(self.data)/math.sqrt(12)
		fit = local_gauss_fit(self.data, mu0, sigma0, A=100, xerrors=xerrs, yerrors=self.errors)
		if plot:
			fit.plot(fit.mu-4*fit.sigma, fit.mu+4*fit.sigma, N=500, color="red")
		return fit

def kalibration(plot=True):
	print("##### CALIBRATION #####")
	channels = []
	channel_errors = []
	energies = []

	with LatexTable("out/calib_peaks.tex") as peaktable, LatexTable("out/calib_counts.tex") as counttable, LatexTable("out/calib_relamp.tex") as relamptable:
		peaktable.header("El.", "Höhe", "Channel", "Breite", r'$\chi^2/\textrm{ndf}$', "Linie", "Energie", lineafter=0)
		peaktable.row("", "", r"$x_\textrm{max}$", r"$\Delta x$", "", "", "")
		peaktable.hline()

		relamptable.header("Element", "Peak", "Relative Höhe", "Literaturwert", lineafter=1)
		counttable.header("Element", "Messdauer", "Ereignisse")

		for filename, meta in kalib.items():
			peaktable.hline()
			relamptable.hline()
			#print("="*10, filename, "="*10)
			experiment = Experiment("data/" + filename, title=meta["title"])
			counttable.row(meta["element"], r'$' + experiment.mcameta["Accumulation Time"] + r' \unit{s}$', experiment.mcameta["Slow Count"])
			if plot:
				plt.clf()
				experiment.errorplot()
				plt.title(meta["title"])
			lines = []
			kalpha = None
			for i, (mu0, sigma0, peak, energy, lrel) in enumerate(sorted(meta["peaks"], key=lambda peak: peak[0]), 1):
				fit = experiment.find_peak(mu0, sigma0, plot=plot)
				channels.append(fit.mu)
				channel_errors.append(fit.sigma_mu)
				#channel_errors.append(fit.sigma) # FEHLER
				energies.append(energy)
				lines.append("Peak {:d}, Channel {:s}, Energy {:s}".format(i, format_error(fit.mu, fit.sigma_mu), format_error(energy, 0.01, unit='keV')))

				if not kalpha:
					kalpha = (fit.A, fit.sigma_A)
				fehler = np.sqrt(np.power(fit.sigma_A/kalpha[0],2)+np.power(fit.A/np.power(kalpha[0],2)*kalpha[1], 2))
				rel = format_error(fit.A / kalpha[0]*100, fehler*100, unit=r'\%')
				if kalpha[0] == fit.A:
					rel = r"100 \%"
				#rel = _n(fit.A/kalpha[0]*100, precision=3) + r"\%"
				lrel = _n(lrel, precision=3) + r"\%"
				relamptable.row((meta["element"], len(meta["peaks"])),
					r'$ ' + peak + r'$', rel, lrel
				)

				peaktable.row((meta["element"], len(meta["peaks"])),
					format_error(fit.A, fit.sigma_A, parenthesis=False),
					format_error(fit.mu, fit.sigma_mu, parenthesis=False),
					format_error(fit.sigma, fit.sigma_sigma, parenthesis=False),
					"%.2f" % fit.chisqndf,
					r'$ ' + peak + r'$',
					r'$ ' + ("%.2f" % energy) + r' \unit{keV}$'
				)
				if plot:
					plt.xlim(fit.mu-5*fit.sigma, fit.mu+5*fit.sigma)
					plt.autoscale(enable=True, axis='y')
					plt.savefig("out/" + clean_filename(filename) + "_peak_%d." % i + SAVETYPE)
			if plot:
				text = "\n".join(lines)
				#info_box(text, location='tr')
				plt.xlim(0,4096)
				plt.autoscale(enable=True, axis='y')
				plt.savefig("out/" + clean_filename(filename) + "_all." + SAVETYPE)

	X = np.array(channels)
	SX = np.array(channel_errors)
	Y = np.array(energies)
	SY = 0.01 # 1 eV

	fit = linear_fit(X, Y, xerr=SX, yerr=SY)
	if plot:
		plt.clf()
		plt.errorbar(X, Y, xerr=SX, yerr=SY, fmt=',')
		fit.plot(np.min(X), np.max(X), box='tl', units={"slope": "keV / Channel", "offset": "keV"})
		plt.xlabel(r"Channel")
		plt.ylabel(r"$ E $ / keV")
		plt.title(r"Kalibrationsgerade: Fit")
		plt.savefig("out/calib_fit." + SAVETYPE)

		plt.clf()
		err = fit.combine_errors(X, SX, SY)
		fit.plot_residuums(X, Y, err, box='bl', fmt=",")
		plt.xlabel(r"Channel")
		plt.ylabel(r"$ E $ / keV")
		plt.title(r"Kalibrationsgerade: Residuen")
		plt.savefig("out/calib_residuum." + SAVETYPE)

	s = format_error(fit.slope*1000, fit.sigma_slope*1000, unit="eV / Kanal", surroundmath=False)
	o = format_error(fit.offset*1000, fit.sigma_offset*1000, unit="eV", surroundmath=False)

	with open('out/calib.tex', 'w') as file:
		file.write(r'Die in diesem Versuch verwendeten Einstellungen f\"uhren zu einer Einteilung von \[ ' +s+ r' \] wobei der erste Kanal die Energie \[ '+o+r' \] besitzt.')

	# print("RESULT: linear fit with chisq: %.2f" % (fit.chisqndf))
	# print("RESULT: energy per bin: %.1f keV" % (fit.slope))
	# print("RESULT: energy of 0th bin: %.1f keV" % (fit.offset))
	return fit


def auswertung(calibration, plot=True, smooth=False):
	print("##### UNKNOWN #####")

	with LatexTable("out/test_peaks.tex") as peaktable, LatexTable("out/test_counts.tex") as counttable:
		peaktable.header("Probe", "Höhe", "Energie", "Breite", lineafter=0)
		peaktable.row("", "", r"$E_\textrm{max}$", r"$\Delta E$")
		peaktable.hline()
		counttable.header("Element", "Messdauer", "Ereignisse", "korrigiert", "Anteil")

		leer = Experiment("data/G20_Leer.mca")
		counttable.row("Leermessung", r'$' + leer.mcameta["Accumulation Time"] + r' \unit{s}$', leer.mcameta["Slow Count"], "-", "-")
		counttable.hline()

		for filename, meta in data.items():
			peaktable.hline()
			#print("="*10, filename, "="*10)
			experiment = Experiment("data/" + filename, title=meta["title"], calibration=calibration)

			count = int(experiment.mcameta["Slow Count"])
			reduced = int(count-0.5*int(leer.mcameta["Slow Count"]))
			percentage = "%.1f \\%%" % (reduced/count * 100)
			counttable.row(meta["title"], r'$' + experiment.mcameta["Accumulation Time"] + r' \unit{s}$', count, reduced, percentage)
			if plot:
				plt.clf()
				experiment.plot()
				plt.title("Probe: " + meta["title"] + " - Rohdaten")
				plt.xlim(0,4096)
				plt.savefig("out/" + clean_filename(filename) + "_raw." + SAVETYPE)

			experiment.subtract_empty("data/G20_Leer.mca", 0.5)

			if smooth:
				experiment.smooth(0.1)

			if plot:
				plt.clf()
				experiment.errorplot()
				plt.title("Probe: " + meta["title"])
				experiment.set_energy_labels(stepsize=0.1)

			lines = []
			for i, (mu0, sigma0) in enumerate(sorted(meta["peaks"], key=lambda peak: peak[0]), 1):
				fit = experiment.find_peak(mu0, sigma0, plot=plot)

				mu = experiment.channel2energy(fit.mu)
				error_mu_stat = experiment.channelwidth2energywidth(fit.sigma_mu)
				error_mu_sys = experiment.sigma_channel2energy(fit.mu)
				string_mu = format_error(mu, error_mu_stat, error_mu_sys, unit='keV')

				sigma = experiment.channelwidth2energywidth(fit.sigma)
				error_sigma_stat =  experiment.channelwidth2energywidth(fit.sigma_sigma)
				#error_sigma_sys =  experiment.sigma_channelwidth2energywidth(fit.sigma) # close to 0
				string_sigma = format_error(sigma, error_sigma_stat, unit='keV')

				lines.append(r"$E_" + str(i) + r"$ = " + string_mu)
				peaktable.row((meta["title"], len(meta["peaks"])),
					format_error(fit.A, fit.sigma_A, parenthesis=False), string_mu, string_sigma
				)

				if plot:
					plt.title("Probe: " + meta["title"] + " - Peak %d" % i)
					plt.xlim(fit.mu-5*fit.sigma, fit.mu+5*fit.sigma)
					plt.autoscale(enable=True, axis='y')
					l, u = plt.ylim()
					plt.ylim(0, u)
					plt.savefig("out/" + clean_filename(filename) + "_peak_%d." % i + SAVETYPE)

			if plot:
				plt.title("Probe: " + meta["title"])
				experiment.set_energy_labels(stepsize=5)
				plt.xlim(0,4096)
				plt.autoscale(enable=True, axis='y')
				text = "\n".join(lines)
				info_box(text, location='tr')
				plt.savefig("out/" + clean_filename(filename) + "_all." + SAVETYPE)


def test_smoothing(calibration):
	print("##### SMOOTHING TEST #####")

	with LatexTable("out/test_smoothing.tex") as peaktable:
		peaktable.header(" ", r'\multicolumn{3}{|c|}{ohne Gl\"attung}', r'\multicolumn{3}{|c|}{mit Gl\"attung}', align=["c"]*7, lineafter=0)
		peaktable.row("Probe", "Höhe", "Breite",  r'$\chi^2/\textrm{ndf}$', "Höhe", "Breite", r'$\chi^2/\textrm{ndf}$')
		peaktable.hline()
		for filename, meta in data.items():
			peaktable.hline(1)

			experiment = Experiment("data/" + filename, calibration=calibration, title=meta["title"])
			experiment2 = Experiment("data/" + filename, calibration=calibration, title=meta["title"])
			experiment.subtract_empty("data/G20_Leer.mca", 0.5)
			experiment2.subtract_empty("data/G20_Leer.mca", 0.5)
			experiment2.smooth(0.1)
			lines = []
			for i, (mu0, sigma0) in enumerate(sorted(meta["peaks"], key=lambda peak: peak[0]), 1):
				fit = experiment.find_peak(mu0, sigma0, plot=plot)
				fit2 = experiment2.find_peak(mu0, sigma0, plot=plot)

				peaktable.row((meta["title"], len(meta["peaks"])),
					_n(fit.mu, precision=3),
					_n(fit.sigma, precision=3),
					_n(fit.chisqndf),
					_n(fit2.mu, precision=3),
					_n(fit2.sigma, precision=3),
					_n(fit2.chisqndf),
				)

def energieaufloesung(calibration, plot=True):
	print("##### ENERGY RESOLUTION #####")

	energies = []
	widths = []
	sigma_energies = []
	sigma_widths = []

	for filename, meta in kalib.items():
		experiment = Experiment("data/" + filename, title=meta["title"], calibration=calibration)
		for peak in meta["peaks"]:
			mu0, sigma0 = peak[0], peak[1]
			fit = experiment.find_peak(mu0, sigma0, plot=False)
			energies.append(experiment.channel2energy(fit.mu))
			widths.append(experiment.channelwidth2energywidth(fit.sigma))
			sigma_energies.append(experiment.channelwidth2energywidth(fit.sigma_mu))
			sigma_widths.append(experiment.channelwidth2energywidth(fit.sigma_sigma))

	for filename, meta in data.items():
		experiment = Experiment("data/" + filename, title=filename, calibration=calibration)
		experiment.subtract_empty("data/G20_Leer.mca", 0.5)
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

	#X  = energies
	#Y = widths/energies
	#SX = sigma_energies
	#SY = Y*np.sqrt(np.power(sigma_widths/widths, 2) + np.power(sigma_energies/energies, 2))
	#func = lambda x, a, b: np.sqrt(np.power(a/np.sqrt(x), 2) + np.power(b,2)) # + np.power(c/x, 2)

	X = energies
	Y = np.power(widths, 2)
	SX = sigma_energies
	SY = 2*widths*sigma_widths

	func = lambda x, a0, a1,: a0+x*a1 #+np.power(x,2)*a2

	fit = Fit(func)
	fit.a0 = 1
	fit.a1 = 1
	#fit.a2 = 1

	for _ in range(10):
		err = fit.combine_errors(X, SX, SY)
		fit.fit(X, Y, err)

	#print("RESULT: error fit with chisq/ndf: %.2f" % (fit.chisqndf))

	if plot:
		plt.clf()
		plt.errorbar(X, Y, xerr=SX, yerr=SY, fmt=',')
		fit.plot(np.min(X), np.max(X), box='br', units={'a0': 'keV', 'a1': r'\sqrt{keV}'})
		plt.xlabel(r"$ E $ / keV")
		plt.ylabel(r"$ \Delta E^2 / keV^2$")
		#plt.title(r'Energieauflösung: Fit zu $ \frac{\Delta E}{E} = \frac{a}{\sqrt{E}} \oplus b $') #  \oplus \frac{c}{E}
		plt.title(r'Energieauflösung: Fit zu $ \Delta E^2 = a_0 + a_1 E $') #  \oplus \frac{c}{E}
		plt.savefig("out/energyresolution_fit." + SAVETYPE)

		plt.clf()
		err = fit.combine_errors(X, SX, SY)
		fit.plot_residuums(X, Y, err, box='tr', fmt=",")
		plt.xlabel(r"$ E $ / keV")
		plt.ylabel(r"$ \Delta E^2 / keV^2$")
		plt.title("Energieauflösung: Residuen")
		plt.savefig("out/energyresolution_residuum." + SAVETYPE)

	with LatexTable("out/energyresolution.tex") as table:
		table.header("Parameter", "Wert")
		table.row("$a_0$", format_error(fit.a0, fit.sigma_a0, unit="keV^2"))
		table.row("$a_1$", format_error(fit.a1, fit.sigma_a1, unit="keV"))
		#table.row("$a_2$", format_error(fit.a2, fit.sigma_a2))
		table.hline()
		a = math.sqrt(fit.a0)
		b = math.sqrt(fit.a1)
		#c = math.sqrt(fit.a2)
		sa = 0.5*fit.sigma_a0/a
		sb = 0.5*fit.sigma_a1/b
		#sc = 0.5*fit.sigma_a2/c
		table.row("$a$", format_error(a, sa, unit="keV"))
		table.row("$b$", format_error(b, sb, unit="\sqrt{keV}"))
		#table.row("$c$", format_error(c, sc))

	with LatexTable("out/energyresolution_examples.tex") as table:
		energies = np.linspace(10,60,6)
		sigmas = np.sqrt(fit.apply(energies))
		deltas = sigmas * 2*math.sqrt(2*math.log(2))
		table.header("Energie", "Auflösung", "Relative Auflösung", "FWHM", lineafter=0)
		table.row("$E$", "$\sigma$", "$\sigma / E$", "$2 \sqrt{2 \ln{2}} \sigma$")
		table.hline(2)
		for energy, sigma, delta in zip(energies, sigmas, deltas):
			e = "%d keV" % energy
			s = _n(sigma*1000) + " eV"
			d = "%d eV" % (delta*1000)
			table.row(e, s, _n(sigma/energy*100) + r"\%", d)

def plot_leermessung(calibration):
	print("##### PLOT EMPTY #####")

	experiment = Experiment("data/G20_Leer.mca", title="Leermessung", calibration=calibration)
	plt.clf()
	experiment.errorplot()
	experiment.set_energy_labels()
	plt.savefig("out/leer." + SAVETYPE)
	#plt.show()

if __name__=="__main__":
	#mng = plt.get_current_fig_manager()
	#mng.window.state('zoomed')

	plot = True

	fit = kalibration(plot)
	plot_leermessung(fit)
	auswertung(fit, plot, False)
	energieaufloesung(fit, plot)
	testsignal()
	test_smoothing(fit)
