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

from funcs import *
from fits import Fit
from number import *

os.chdir(os.path.dirname(os.path.realpath(__file__)))

SAVETYPE = "png"

def sort_multiple(keys, *others):
	order = keys.argsort()
	result = [keys[order]]
	for other in others:
		result.append(other[order])
	return tuple(result)

def plot_temperature():
	temperature, voltage, error_voltage = np.loadtxt("data/temperatur.txt", unpack=True)
	error_temperature = 0.1 / math.sqrt(12)

	plt.clf()
	plt.errorbar(temperature, voltage, xerr=error_temperature, yerr=error_voltage, fmt='-')
	plt.xlabel("Temperatur / Grad")
	plt.ylabel("Spannung / V")
	plt.title("Temperaturabhängigkeit des Diodenlasers")
	plt.show()

def plot_diode_kennlinie():
	current, power = np.loadtxt("data/kennlinie_diode.txt", unpack=True)
	error_current = 1 / math.sqrt(12)
	error_power = 1 / math.sqrt(12)

	fit_indices_list = []
	fit_indices_list.append(power > 10)
	fit_indices_list.append(np.logical_and(power > 10, current < 420))
	fit_indices_list.append(np.logical_and(power > 10, current > 420))

	fit = Fit(LINEAR)
	fit.set_params(slope=0.85, offset=-160)
	fit.set_data(xdata=current, ydata=power, xerrors=error_current, yerrors=error_power)
	fit.set_labels(xlabel="Strom / mA", ylabel="Leistung / mW")

	for fit_indices in fit_indices_list:
		subfit = fit.filtered(fit_indices)

		plt.clf()
		plt.errorbar(current, power, xerr=error_current, yerr=error_power, fmt=',', color="black")
		subfit.iterative_fit(5)
		subfit.plot(box="tl", color="red", plot_data=False)
		plt.show()

		plt.clf()
		subfit.plot_residual(box="br", color="black", fmt="s")
		plt.show()

		i_threshold = -subfit.ufloat("offset")/subfit.ufloat("slope")
		print("Threshold:", formatUFloat(i_threshold, unit="mA"))

def plot_yag_kennlinie():
	current, power, error_power = np.loadtxt("data/kennlinie_yag2.txt", unpack=True)
	error_current = 1 / math.sqrt(12)

	fit_indices_list = []
	fit_indices_list.append(power > 0.10)
	fit_indices_list.append(np.logical_and(power > 0.10, current < 480))
	fit_indices_list.append(np.logical_and(power > 0.10, current > 480))

	fit = Fit(POLY2)
	fit.set_data(xdata=current, ydata=power, xerrors=error_current, yerrors=error_power)
	fit.set_labels(xlabel="Strom / mA", ylabel="Leistung / mW")

	for fit_indices in fit_indices_list:
		plt.clf()
		fit.plot(plot_fit=False)
		subfit = fit[fit_indices]
		subfit.iterative_fit(5)
		subfit.plot(plot_data=False, box="tl")
		plt.show()

		plt.clf()
		subfit.plot_residual(box="br", color="black", fmt="s")
		plt.show()

def plot_qswitch():
	frequency, error_frequency, power = np.loadtxt("data/frequenz_qswitch.txt", unpack=True)
	error_power = 1
	plt.errorbar(frequency, power, xerr=error_frequency, yerr=error_power, fmt=',')
	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Leistung / uW")
	plt.show()

if __name__=="__main__":
	#plot_temperature()
	plot_diode_kennlinie()
	#plot_yag_kennlinie()
	#plot_qswitch()
