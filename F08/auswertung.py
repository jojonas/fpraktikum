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

def current2diode_power(current, error_current):
	lower_threshold = 182.10
	error_lower_threshold = 0.23
	upper_threshold = 201.6
	error_upper_threshold = 0.5
	lower_slope = 0.8139
	error_lower_slope = 0.0013
	upper_slope = 0.8892
	error_upper_slope = 0.0014

	lower = np.logical_and(current > lower_threshold, current < 410)
	upper = current >= 410

	diode_power = np.zeros_like(current)
	diode_power[lower] = (current[lower] - lower_threshold) * lower_slope
	diode_power[upper] = (current[upper] - upper_threshold) * upper_slope

	error_diode_power = np.zeros_like(diode_power)
	error_diode_power[lower] = np.sqrt(np.power(error_current*lower_slope,2) + np.power((current[lower]-lower_threshold)*error_lower_slope,2) + np.power(lower_slope*error_lower_threshold,2))
	error_diode_power[upper] = np.sqrt(np.power(error_current*upper_slope,2) + np.power((current[upper]-upper_threshold)*error_upper_slope,2) + np.power(upper_slope*error_upper_threshold,2))

	return diode_power, error_diode_power

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
	#fit_indices_list.append(np.logical_and(power > 0.10, lower))
	#fit_indices_list.append(np.logical_and(power > 0.10, upper))

	diode_power, error_diode_power = current2diode_power(current, error_current)

	zero_power = power[current < 150]
	background = zero_power.mean()
	print("Background:", background, "mW")
	power -= background

	plt.clf()
	plt.errorbar(diode_power, power, xerr=error_diode_power, yerr=error_power, fmt=',', color="black")
	plt.xlabel("Diodenleistung / mW")
	plt.ylabel("Laserleistung / mW")
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


		print(subfit)

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

def diode_power2yag_power(diode_power, error_diode_power):
	a0 = -0.272
	error_a0 = 0.005
	a1 = 0.01605
	error_a1 = 0.00010
	a2 = 1.69e-5
	error_a2 = 0.04e-5

	threshold = 16.68

	yag_power = np.power(diode_power, 2)*a2 + diode_power*a1 + a0
	yag_power[diode_power < threshold] = 0
	error_yag_power = np.sqrt(np.power((2*diode_power*a2+a1)*error_diode_power, 2) + np.power(np.power(diode_power,2)*error_a2, 2) + np.power(diode_power*error_a1, 2) + np.power(error_a0, 2))
	error_yag_power[diode_power < threshold] = 0

	return yag_power, error_yag_power

def plot_ktp():
	current, power, error_power = np.loadtxt("data/ktp_kristall.txt", unpack=True)
	error_current = 1 / math.sqrt(12)

	# uW -> mW
	power /= 1000
	error_power /= 1000

	plt.clf()
	plt.errorbar(current, power, xerr=error_current, yerr=error_power, fmt=',', color="black")
	plt.xlabel("Strom / mA")
	plt.ylabel("Leistung mit KTP / mW")
	plt.xlim(0,700)
	plt.savefig("out/ktp_raw." + SAVETYPE)

	lower = np.logical_and(current > 182.10, current < 410)
	upper = current >= 410

	diode_power, error_diode_power = current2diode_power(current, error_current)

	plt.clf()
	plt.errorbar(diode_power, power, xerr=error_diode_power, yerr=error_power, fmt=',', color="black")
	plt.xlabel("Diodenleistung / mW")
	plt.ylabel("Laserleistung mit KTP / mW")
	plt.savefig("out/ktp_raw2." + SAVETYPE)

	yag_power, error_yag_power = diode_power2yag_power(diode_power, error_diode_power)

	plt.clf()
	plt.errorbar(yag_power, power, xerr=error_yag_power, yerr=error_power, fmt=',', color="black")
	plt.xlabel("Laserleistung ohne KTP / mW")
	plt.ylabel("Laserleistung mit KTP / mW")
	plt.xlim(0,10)
	plt.savefig("out/ktp_raw3." + SAVETYPE)

	plt.clf()
	fit = Fit(POLY2)
	fit.set_data(xdata=yag_power, ydata=power, xerrors=error_yag_power, yerrors=error_power)
	fit.set_labels(xlabel="Laserleistung ohne KTP / mW", ylabel="Laserleistung mit KTP / mW")

	fit = fit.filtered(yag_power > 2)

	fit.iterative_fit(5)

	fit.plot(box="tl", units={"a2": "1/mW", "a0": "mW"})
	plt.savefig("out/ktp_fit." + SAVETYPE)

	plt.clf()
	fit.plot_residual(box="tl")
	plt.savefig("out/ktp_residual." + SAVETYPE)

def plot_qswitch():
	frequency, error_frequency, power = np.loadtxt("data/frequenz_qswitch.txt", unpack=True)
	error_power = 1 / math.sqrt(12)

	power /= 1000
	error_power /= 1000

	plt.clf()
	plt.errorbar(frequency, power, xerr=error_frequency, yerr=error_power, fmt=',', color="black")
	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Leistung / mW")
	plt.savefig("out/qswitch_raw." + SAVETYPE)

	#power -= power[0]

	fit = Fit(POLY3)
	fit.set_data(xdata=frequency, ydata=power, xerrors=error_frequency, yerrors=error_power)
	fit.set_params(a0=380, a1=0.005, a2=-7e-8)
	fit.set_labels(xlabel="Frequenz / Hz", ylabel="Leistung / mW")

	subfit = fit.filtered(np.logical_and(frequency < 20e3, frequency > 1e3))
	subfit.iterative_fit(5)

	plt.clf()
	subfit.plot(box="br", units={"a3": "mW/Hz^3", "a2": "mW/Hz^2", "a1": "mW/Hz", "a0": "mW"})
	plt.savefig("out/qswitch_power_fit." + SAVETYPE)

	plt.clf()
	subfit.plot_residual(box="tr")
	plt.savefig("out/qswitch_power_residual." + SAVETYPE)

	plt.clf()
	time, voltage = _load_oscilloscope_csv("data/ALL0023/F0023CH2.CSV")
	time *= 1E6
	time -= -935

	plt.plot(time, voltage, ".", color="black")

	pause = np.logical_and(time > 300, time < 900)
	one_period = np.logical_and(time > 0, time < 1000)
	print(one_period, one_period.sum())

	mean_voltage = ufloat(voltage[one_period].mean(), voltage[pause].std() / math.sqrt(one_period.sum()))
	peak_voltage = ufloat(voltage[one_period].max(), voltage[pause].std())

	print("ErrorL:", voltage[pause].std())

	print("Mean voltage:", mean_voltage, "V")
	print("Peak voltage:", peak_voltage, "V")

	frequency = 1.01798e3 # Hz
	mean_power = subfit.ueval(frequency)
	print("Mean power:", mean_power*1000, "uW")

	peak_power = mean_power * peak_voltage / mean_voltage
	print("Peak power:", peak_power*1000, "uW")

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
	#plot_spiking()
	#plot_ktp()
	plot_qswitch()
