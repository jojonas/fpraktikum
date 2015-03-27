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
from uncertainties import ufloat

from funcs import *
from fits import Fit, info_box
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

	plt.clf()
	plt.errorbar(current, power, xerr=error_current, yerr=error_power, fmt=',', color="black")
	plt.xlabel("Strom / mA")
	plt.ylabel("Leistung / mW")
	plt.xlim(0,700)
	plt.savefig("out/kennlinie_yag_raw." + SAVETYPE)

	lower = np.logical_and(current > 182.10, current < 410)
	upper = current >= 410

	fit_indices_list = []
	fit_indices_list.append(power > 0.10)
	fit_indices_list.append(np.logical_and(power > 0.10, lower))
	fit_indices_list.append(np.logical_and(power > 0.10, upper))

	diode_power = np.zeros_like(current)
	diode_power[lower] = (current[lower] - 182.10) * 0.8139
	diode_power[upper] = (current[upper] - 201.6) * 0.8892

	error_diode_power = np.zeros_like(diode_power)
	error_diode_power[lower] = np.sqrt(np.power(error_current*0.8139,2) + np.power((current[lower]-182.10)*0.0013,2) + np.power(0.8139*0.23,2))
	error_diode_power[upper] = np.sqrt(np.power(error_current*0.8892,2) + np.power((current[upper]-201.6)*0.0014,2) + np.power(0.8892*0.5,2))

	zero_power = power[current < 150]
	background = zero_power.mean()
	print("Background:", background, "mW")
	power -= background


	plt.clf()
	plt.errorbar(diode_power, power, xerr=error_diode_power, yerr=error_power, fmt=',', color="black")
	plt.xlabel("Diodenleistung / mW")
	plt.ylabel("laserleistung / mW")
	plt.savefig("out/kennlinie_yag_raw2." + SAVETYPE)

	func = lambda x, slope, threshold: (x - threshold) * slope
	fit = Fit(func)
	fit.set_params(slope=0.016, threshold=16)
	fit.set_data(xdata=diode_power, ydata=power, xerrors=error_diode_power, yerrors=error_power)
	fit.set_labels(xlabel="Diodenleistung / mW", ylabel="Laserleistung / mW")
	fit.iterative_fit(1)

	plt.clf()
	fit.plot(plot_data=True, plot_fit=True, box="tr")
	plt.savefig("out/kennlinie_yag_linear_fit." + SAVETYPE)
	plt.clf()
	fit.plot_residual(box="br", fmt="s")
	plt.savefig("out/kennlinie_yag_linear_residual." + SAVETYPE)

	fit = Fit(POLY2)
	fit.set_data(xdata=diode_power, ydata=power, xerrors=error_diode_power, yerrors=error_power)
	fit.set_labels(xlabel="Diodenleistung / mW", ylabel="Laserleistung / mW")
	fit.iterative_fit(5)

	for i, fit_indices in enumerate(fit_indices_list):
		plt.clf()
		fit.plot(plot_fit=False)
		subfit = fit[fit_indices]
		subfit.iterative_fit(5)

		#zero_power = power[power < 0.03]
		#zero = umean(zero_power)
		#plt.axhline(zero.n, color="blue")

		x = 200
		try:
			i_threshold = solve_quadratic(subfit.uvalue("a0"), subfit.uvalue("a1"), subfit.uvalue("a2"))
			i_slope = 2*x*subfit.uvalue("a2") + subfit.uvalue("a1")
		except ValueError:
			print("no solution", i)
		else:
			print("Threshold:", formatUFloat(i_threshold, unit="mW"))
			print("Efficiency:", formatUFloat(i_slope*100, unit="%"))
			lines = "threshold = " + formatUFloat(i_threshold, unit="mW") + "\n" + "efficiency at %d mW = " % x + formatUFloat(i_slope*100, unit="%")
			info_box(lines, location="tl")

		subfit.plot(plot_data=False)
		plt.savefig("out/kennlinie_yag_%d_fit." % i + SAVETYPE)

		plt.clf()
		subfit.plot_residual(box="br", color="black", fmt="s")
		plt.savefig("out/kennlinie_yag_%d_residual." % i + SAVETYPE)

		plt.clf()
		x = diode_power[diode_power > i_threshold.n]
		plt.plot(x, 100*(2*x*subfit.value("a2") + subfit.value("a1")), color="black")
		#x_new = diode_power - i_threshold.n
		#plt.plot(diode_power, 100*(x_new*subfit.value("a2") + subfit.value("a1") + subfit.value("a0")/x_new), color="red")
		plt.xlabel("Diodenleistung / mW")
		plt.ylabel("Effizienz / %")
		plt.savefig("out/kennlinie_yag_%d_efficiency." % i + SAVETYPE)

def plot_spiking():
	time, voltage = _load_oscilloscope_csv("data/ALL0015/F0015CH2.CSV")
	time -= 0.003297480000
	time *= 1E6

	voltage -= voltage[time < 0].mean()

	plt.clf()
	plt.plot(time, voltage, ".", color="black")
	plt.xlabel("Zeit / us")
	plt.ylabel("Spannung / V")
	plt.xlim(-10, 80)
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
	#plot_yag_lifetime()
	#plot_yag_kennlinie()
	plot_spiking()
	#plot_qswitch()
