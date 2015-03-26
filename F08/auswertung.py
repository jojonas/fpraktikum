import os
import math

from pylab import rcParams
rcParams['savefig.dpi'] = 200
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

from uncertainties import wrap as uwrap
from uncertainties import ufloat

from funcs import *
from fits import Fit
from number import *

os.chdir(os.path.dirname(os.path.realpath(__file__)))

SAVETYPE = "png"

def _load_oscilloscope_csv(filename):
	return np.loadtxt(filename, skiprows=18, usecols=(3,4), delimiter=",", unpack=True)

def sort_multiple(keys, *others):
	order = keys.argsort()
	result = [keys[order]]
	for other in others:
		result.append(other[order])
	return tuple(result)

@uwrap
def solve_quadratic(a0, a1, a2):
	p = a1/a2
	q = a0/a2
	a=-p/2
	b=math.sqrt(math.pow(p/2, 2)-q)
	return a+b

def umean(data):
	N = len(data)
	mean = data.mean()
	std = data.std()
	return ufloat(mean, std/math.sqrt(N))

def plot_temperature():
	temperature, voltage, error_voltage = np.loadtxt("data/temperatur.txt", unpack=True)
	error_temperature = 0.1 / math.sqrt(12)

	plt.clf()
	plt.errorbar(temperature, voltage, xerr=error_temperature, yerr=error_voltage, fmt=',', color="black")
	# plt.axvline(25.5, color="gray")
	# plt.axvline(29.0, color="gray")
	# plt.axvline(32.3, color="gray")
	# plt.axvline(34.8, color="gray")
	# plt.axvline(37.5, color="gray")
	# plt.axvline(40.6, color="gray")
	plt.xlabel("Temperatur / Grad")
	plt.ylabel("Spannung / V")
	plt.savefig("out/temperatur." + SAVETYPE)


def plot_diode_kennlinie():
	current, power = np.loadtxt("data/kennlinie_diode.txt", unpack=True)
	error_current = 1 / math.sqrt(12)
	error_power = 1 / math.sqrt(12)

	plt.clf()
	plt.errorbar(current, power, xerr=error_current, yerr=error_power, fmt=',', color="black")
	plt.xlabel("Strom / mA")
	plt.ylabel("Leistung / mW")
	plt.savefig("out/kennlinie_diode_raw." + SAVETYPE)

	fit_indices_list = []
	fit_indices_list.append(current > 200)

	valid = current > 200
	valid[current > 410] = False
	fit_indices_list.append(valid)

	valid = current > 410
	fit_indices_list.append(valid)

	zero_power = power[current < 150]
	background = zero_power.mean()
	print("Background:", background, "mW")
	power -= background

	func = lambda x, slope, threshold: (x - threshold) * slope
	fit = Fit(func)
	fit.set_params(slope=0.85, threshold=190)
	fit.set_data(xdata=current, ydata=power, xerrors=error_current, yerrors=error_power)
	fit.set_labels(xlabel="Strom / mA", ylabel="Leistung / mW")

	for i, fit_indices in enumerate(fit_indices_list):
		subfit = fit.filtered(fit_indices)

		plt.clf()
		plt.errorbar(current, power, xerr=error_current, yerr=error_power, fmt=',', color="black")
		subfit.iterative_fit(5)

		#plt.axhline(background, color="blue")
		#plt.axvline(subfit.value("threshold"), color="blue")
		print("Threshold:", formatUFloat(subfit.uvalue("threshold"), unit="mA"))
		print("Slope:", formatUFloat(subfit.uvalue("slope"), unit="W/A"))
		print("Chi^2/ndf:", subfit.chi2ndf)
		print("===")

		subfit.plot(range=(subfit.value("threshold"), current.max()), color="black", plot_fit=False, plot_data=True, fmt="s")
		subfit.plot(range=(subfit.value("threshold"), current.max()), box="tl", color="red", plot_data=False, units={"threshold": "mA", "slope": "W/A"})
		plt.savefig("out/kennlinie_diode_%d_fit." % i + SAVETYPE)

		plt.clf()
		subfit.plot_residual(box="br", color="black", fmt="s")
		plt.savefig("out/kennlinie_diode_%d_residual." % i + SAVETYPE)


def plot_yag_lifetime():
	time, voltage = _load_oscilloscope_csv("data/ALL0010/F0010CH2.CSV")
	time *= 1E6 # us
	voltage *= 1000 # mV

	error_voltage = voltage[time<0].std()
	print("Fehler:", error_voltage, "mV")

	fit = Fit(EXPONENTIAL_DECAY)
	fit.set_data(xdata=time, ydata=voltage, yerrors=error_voltage)
	fit.set_labels(xlabel="Zeit / us", ylabel="Spannung / mV")
	fit = fit.filtered(np.logical_and(time>0, time<1500))
	fit.set_params(A=0.1, offset=0.20, T=250)
	fit.iterative_fit(1)

	print("Lebensdauer:", formatUFloat(fit.uvalue("T"), "us"), "- Chi^2:", fit.chi2, "- ndf:", fit.ndf, "- Chi^2/ndf:", fit.chi2ndf)

	plt.clf()
	fit.plot(box="tr", units={"T": "us", "A": "mV", "offset": "mV"})
	plt.savefig("out/yag_lifetime." + SAVETYPE)

def plot_yag_kennlinie():
	current, power, error_power = np.loadtxt("data/kennlinie_yag2.txt", unpack=True)
	error_current = 1 / math.sqrt(12)

	fit_indices_list = []
	fit_indices_list.append(power > 0.10)
	fit_indices_list.append(np.logical_and(power > 0.10, current < 480))

	fit = Fit(POLY2)
	fit.set_data(xdata=current, ydata=power, xerrors=error_current, yerrors=error_power)
	fit.set_labels(xlabel="Strom / mA", ylabel="Leistung / mW")

	for i, fit_indices in enumerate(fit_indices_list):
		plt.clf()
		fit.plot(plot_fit=False)
		subfit = fit[fit_indices]
		subfit.iterative_fit(5)

		zero_power = power[power < 0.03]
		zero = umean(zero_power)
		plt.axhline(zero.n, color="blue")

		i_threshold = solve_quadratic(subfit.uvalue("a0")-zero, subfit.uvalue("a1"), subfit.uvalue("a2"))
		plt.axvline(i_threshold.n, color="blue")
		print("Threshold:", formatUFloat(i_threshold, unit="mA"))

		subfit.plot(plot_data=False, box="tl")
		plt.savefig("out/kennlinie_yag_%d_fit." % i + SAVETYPE)

		plt.clf()
		subfit.plot_residual(box="br", color="black", fmt="s")
		plt.savefig("out/kennlinie_yag_%d_residual." % i + SAVETYPE)


def plot_spiking():
	time, voltage = _load_oscilloscope_csv("data/ALL0015/F0015CH2.CSV")
	time *= 1E6
	time -= time[0]

	plt.clf()
	plt.plot(time, voltage, ".", color="black")
	plt.xlabel("Zeit / us")
	plt.ylabel("Spannung / V")
	plt.savefig("out/spiking." + SAVETYPE)

def plot_qswitch():
	frequency, error_frequency, power = np.loadtxt("data/frequenz_qswitch.txt", unpack=True)
	error_power = 1

	fit = Fit(POLY2)
	fit.set_data(xdata=frequency, ydata=power, xerrors=error_frequency, yerrors=error_power)
	fit.set_params(a0=380, a1=0.005, a2=-7e-8)
	fit.set_labels(xlabel="Frequenz / Hz", ylabel="Leistung / uW")

	subfit = fit.filtered(np.logical_and(frequency < 20e3, frequency > 1e3))
	subfit.iterative_fit(5)

	plt.clf()
	subfit.plot() #box="br")
	plt.savefig("out/qswitch_power_fit." + SAVETYPE)

	plt.clf()
	subfit.plot_residual(box="tr")
	plt.savefig("out/qswitch_power_residual." + SAVETYPE)

	plt.clf()
	time, voltage = _load_oscilloscope_csv("data/ALL0023/F0023CH2.CSV")
	time *= 1E6
	time -= time[0]

	plt.plot(time, voltage, ".", color="black")

	voltage_noise = voltage[np.logical_and(time > 500, time < 1000)].std()

	mean_voltage = ufloat(voltage.mean(), voltage_noise)
	peak_voltage = ufloat(voltage.max(), voltage_noise)

	frequency = 1.01798e3 # Hz
	mean_power = subfit.ueval(frequency)

	peak_power = mean_power * peak_voltage / mean_voltage
	print("Peak power:", formatUFloat(peak_power), "uW")

	plt.axhline(mean_voltage.n, color="red", linewidth=2)
	plt.axhline(peak_voltage.n, color="red", linewidth=2)

	plt.xlabel("Zeit / us")
	plt.ylabel("Spannung / V")

	plt.savefig("out/qswitch_osci." + SAVETYPE)


if __name__=="__main__":
	#plot_temperature()
	#plot_diode_kennlinie()
	plot_yag_lifetime()
	#plot_yag_kennlinie()
	#plot_spiking()
	#plot_qswitch()
