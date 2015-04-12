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

from scipy.interpolate import interp1d

os.chdir(os.path.dirname(os.path.realpath(__file__)))

SAVETYPE = "png"

def siPrefix(number):
	s = str(number)
	table = {
		r'e-12': 'f',
		r'e-09': 'n',
		r'e-06': 'u',
		r'e-03': 'm',
		r'e+03': 'k',
		r'e+06': 'M',
		r'e+09': 'G',
		r'e+12': 'T',
		r'e+15': 'P',
	}
	for suffix, prefix in table.items():
		if s.endswith(suffix):
			s = s[:-len(suffix)] + " " + prefix
	return s

def uvline(x, **opts):
	plt.axvline(x.n, linewidth=2, **opts)
	plt.axvline(x.n-x.s, linewidth=1, **opts)
	plt.axvline(x.n+x.s, linewidth=1, **opts)

def uhline(x, **opts):
	plt.axhline(x.n, linewidth=2, **opts)
	plt.axhline(x.n-x.s, linewidth=1, **opts)
	plt.axhline(x.n+x.s, linewidth=1, **opts)

def plot_teil1():
	C = 94.5e-9 # F
	R = 1.487e3 # Ohm

	freq, multimeter_voltage, oscilloscope_voltage = np.loadtxt("data/teil1b.txt", unpack=True)

	plt.clf()
	plt.minorticks_on()
	plt.plot(freq, multimeter_voltage, "s", color="gray", label="Multimeter")
	plt.plot(freq, oscilloscope_voltage, "o", color="red", label="Oszilloskop")
	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Spannung / V")
	plt.xscale("log")
	plt.legend()
	plt.savefig("out/tiefpass_raw." + SAVETYPE)

	print("BAR MULTI", multimeter_voltage[:10].mean())
	print("BAR OSZI", oscilloscope_voltage[:10].mean() / math.sqrt(2))

	multimeter_voltage /= multimeter_voltage[:10].mean()
	oscilloscope_voltage /= oscilloscope_voltage[:10].mean()

	f_grenz = 1 / (2*math.pi * R * C)

	plt.clf()
	plt.minorticks_on()
	plt.plot(freq, multimeter_voltage, "s", color="gray", label="Multimeter")
	plt.plot(freq, oscilloscope_voltage, "o", color="red", label="Oszilloskop")
	plt.axvline(f_grenz, color="black", linestyle="--")
	plt.axhline(1 / math.sqrt(2), color="black", linestyle="--")

	theorie = 1.000 / np.sqrt(1 + np.power(2*np.pi*freq*C*R, 2))
	plt.plot(freq, theorie, "-", color="black", label="Theorie")

	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Ua / Ue")
	plt.xscale("log")
	plt.legend()
	plt.savefig("out/tiefpass." + SAVETYPE)

def plot_tiefpass():
	freq, error_freq, multimeter_voltage, oscilloscope_voltage = np.loadtxt("data/teil2_tiefpass.txt", unpack=True)

	plt.clf()
	plt.minorticks_on()
	plt.plot(freq, multimeter_voltage, "s", color="gray", label="Multimeter")
	plt.plot(freq, oscilloscope_voltage, "o", color="red", label="Oszilloskop")
	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Spannung / V")
	plt.xscale("log")
	plt.legend()
	plt.savefig("out/li_tiefpass_raw." + SAVETYPE)

	full = oscilloscope_voltage[:2].mean()
	#oscilloscope_voltage /= full

	threshold = full / math.sqrt(2)

	f_grenz = ufloat(642, 30)

	plt.clf()
	plt.minorticks_on()
	#plt.plot(freq, multimeter_voltage, "s")
	plt.plot(freq, oscilloscope_voltage, "o", color="red")

	plt.axhline(threshold, color="black", linestyle="--")
	plt.axvline(f_grenz.n, color="red")
	plt.axvline(1000, color="black", linestyle="--")
	plt.axhline(0, color="black", linestyle="--")
	info_box("f0 = 1000 Hz\nfg = ({:P}) Hz".format(f_grenz), location="tr")

	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Spannung / V")
	plt.xscale("log")

	plt.savefig("out/li_tiefpass." + SAVETYPE)

def plot_hochpass():
	freq, error_freq, multimeter_voltage, oscilloscope_voltage = np.loadtxt("data/teil2_hochpass.txt", unpack=True)

	plt.clf()
	plt.minorticks_on()
	plt.plot(freq, multimeter_voltage, "s", color="gray", label="Multimeter")
	plt.plot(freq, oscilloscope_voltage, "o", color="red", label="Oszilloskop")
	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Spannung / V")
	plt.xscale("log")
	plt.legend()
	plt.savefig("out/li_hochpass_raw." + SAVETYPE)

	full = oscilloscope_voltage[-2:].mean()
	#oscilloscope_voltage /= full

	threshold = full / math.sqrt(2)

	f_grenz = ufloat(1273, 20)

	plt.clf()
	plt.minorticks_on()
	#plt.plot(freq, multimeter_voltage, "s")
	plt.plot(freq, oscilloscope_voltage, "o", color="red")

	plt.axhline(threshold, color="black", linestyle="--")
	plt.axvline(f_grenz.n, color="red")
	plt.axvline(1000, color="black", linestyle="--")
	plt.axhline(0, color="black", linestyle="--")
	info_box("f0 = 1000 Hz\nfg = ({:P}) Hz".format(f_grenz), location="tr")

	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Spannung / V")
	plt.xscale("log")
	plt.savefig("out/li_hochpass." + SAVETYPE)

def plot_bandpass():
	freq, error_freq, multimeter_voltage, oscilloscope_voltage = np.loadtxt("data/teil2_bandpass.txt", unpack=True)
	plt.clf()
	plt.minorticks_on()
	plt.plot(freq, multimeter_voltage, "s", color="gray", label="Multimeter")
	plt.plot(freq, oscilloscope_voltage, "o", color="red", label="Oszilloskop")
	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Spannung / V")
	plt.xscale("log")
	plt.legend()
	plt.savefig("out/li_bandpass_raw." + SAVETYPE)

	#plt.plot(freq, multimeter_voltage, "s")
	#plt.plot(freq, oscilloscope_voltage, "o")

	#func = lambda x, s, t, sigma, l, g: l*LORENTZ(np.log10(x), s, t) + g*GAUSS(np.log10(x), t, sigma, 1)
	func = lambda x, s, t: LORENTZ(np.log10(x), s, t)

	fit = Fit(func)
	fit.set_data(xdata=freq, xerrors=error_freq, ydata=oscilloscope_voltage, yerrors=0.01)
	fit.set_labels(xlabel="Frequenz / Hz", ylabel="Spannung / V")
	#fit.set_params(s=1, t=3, sigma=3, l=1, g=0.5)
	fit.set_params(s=1, t=1000)
	fit.iterative_fit(5)
	plt.clf()
	plt.minorticks_on()
	fit.plot(box="tl", logx=True)
	plt.xlim(100, 10000)
	plt.savefig("out/li_bandpass_fit." + SAVETYPE)
	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tl", logx=True)
	plt.xlim(100, 10000)
	plt.savefig("out/li_bandpass_residual." + SAVETYPE)

	s = fit.uvalue("s")
	t = fit.uvalue("t")

	full = 1/(math.pi * s)
	f0 = umath.pow(10, t)
	threshold = full / math.sqrt(2)

	a = t
	b = math.sqrt(math.sqrt(2) - 1) * s
	f_grenz_lower = umath.pow(10, a-b)
	f_grenz_upper = umath.pow(10, a+b)

	#f_grenz_lower = ufloat(708, 3)
	#f_grenz_upper = ufloat(1390, 5)
	#f_grenz_lower = ufloat(590, 10)
	#f_grenz_upper = ufloat(1665, 25)

	B = (f_grenz_upper - f_grenz_lower)
	Q = f0 / B

	plt.clf()
	plt.minorticks_on()
	#plt.plot(freq, multimeter_voltage, "s")
	plt.plot(freq, oscilloscope_voltage, "o", color="red")

	uhline(full, color="black", linestyle="--")
	uhline(threshold, color="black", linestyle="--")
	uvline(f_grenz_lower, color="red")
	uvline(f_grenz_upper, color="red")
	plt.axvline(1000, color="red", linestyle="-")
	plt.axhline(0, color="black", linestyle="--")
	info_box("f0 = ({:P}) Hz\nf1 = ({:P}) Hz\nf2 = ({:P}) Hz\n\nU(f0) = ({:P}) V\nB = ({:P}) Hz\nQ = {:P}".format(f0, f_grenz_lower, f_grenz_upper, full, B, Q), location="bl")

	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Spannung / V")
	plt.xscale("log")
	#plt.show()
	plt.savefig("out/li_bandpass." + SAVETYPE)


def plot_integration_time():
	freq, mean, pk2pk = np.loadtxt("data/teil4.txt", unpack=True)
	mean = -mean

	mean /= 50 # Verstaerkung
	mean *= 1000 # V->mV
	pk2pk /= 50 # Verstaerkung
	pk2pk *= 1000 # V->mV

	plt.clf()
	plt.minorticks_on()
	plt.errorbar(freq, mean, yerr=pk2pk, fmt="s")
	plt.xlabel("Frequenz / Hz")
	plt.ylabel("Spannung / mV")
	plt.savefig("out/integration_time." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	plt.errorbar(freq/100, np.log(mean), yerr=pk2pk/mean, fmt="s")
	fit = Fit(LINEAR)
	fit.set_data(freq/100, np.log(mean), yerrors=pk2pk/mean)
	fit.set_labels(xlabel="Perioden", ylabel="ln(Spannung / mV)")
	fit = fit.filtered(np.logical_and(freq <= 700, freq >= 100))
	fit.iterative_fit(5)
	fit.plot(box="bl", color="red")
	plt.savefig("out/integration_time_log_fit." + SAVETYPE)
	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tl")
	plt.savefig("out/integration_time_log_residual." + SAVETYPE)

	print("Vorfaktor:", umath.exp(fit.uvalue("offset")))

def faktor():
	fR = ufloat(980, 3)
	I0 = 1e-3 # 1 mA

	durchmesser = ufloat(2.9e-3, 0.05e-3)
	lange_seite = ufloat(3e-3, 0.05e-3)
	kurze_seite = ufloat(2.55e-3, 0.05e-3)
	Vprobe = math.pi * (durchmesser/2)**2 * (lange_seite + kurze_seite)/2

	mu0 = 4*math.pi * 1e-7

	Nsek = 700
	Nprim = 2046
	lsek = 50 / 1000 # 50mm
	lprim = 20 / 1000 # 20mm

	C = 2*math.pi*fR*mu0*I0*Nsek*Nprim/(lsek*lprim)

	chi1 = 5
	nM = ufloat(0.27, 0.01)
	s  = C*chi1*Vprobe*(1-nM)

	v = 10/s

	#S = "{:g}".format(s)
	#V = "{:S}".format(v)

	print("Minimale Empfindlichkeit:", s, "V")
	print("Maximale Verstärkung: ", v)

def load_samples(filename):
	comma_fix = lambda s: float(s.decode("utf-8").strip().replace(",", "."))
	return np.loadtxt(filename, converters={0: comma_fix, 1: comma_fix}, unpack=True)

def plot_cold_empty(error_temperature, error_voltage):
	temperature, voltage = load_samples("data/aufwaermen_ohne.dat")
	time = np.arange(0, len(temperature)) / 3

	plt.clf()
	plt.minorticks_on()
	fit = Fit(POLY4)
	fit.set_data(xdata=temperature, xerrors=error_temperature, ydata=voltage, yerrors=error_voltage)
	fit.set_labels("Temperatur / K", "gemessene Spannung / V")
	fit = fit.filtered(np.logical_and(time > 30, temperature < 120))
	fit.iterative_fit(5)
	fit.plot(box=None)
	twin = secondary_y_axis(lambda x: x * 100)
	twin.set_ylabel("Signalspannung / uV")
	plt.savefig("out/cold_empty_fit." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tl")
	twin = secondary_y_axis(lambda x: x * 100)
	twin.set_ylabel("Signalspannung / uV")
	plt.savefig("out/cold_empty_residual." + SAVETYPE)

	return fit

def plot_cold_noise():
	# Realteil Chi'
	temperature, voltage = load_samples("data/abkuehlen_supra.dat")
	time = np.arange(0, len(temperature)) / 3

	T = 1000

	plt.clf()
	plt.minorticks_on()
	plt.plot(time[time>T], temperature[time>T], ".")
	plt.xlabel("Zeit / s")
	plt.ylabel("Temperatur / K")

	fit = Fit(LINEAR)
	fit.set_data(time, temperature, yerrors=0.01)
	fit = fit.filtered(time > T)
	fit.iterative_fit(1)
	fit.plot(plot_data=False)
	plt.xlim(T,time.max())
	plt.savefig("out/temperature_noise." + SAVETYPE)

	temperature = temperature - fit.eval(time)
	error_temperature = temperature[time > T].std()
	print("Temperature noise:", error_temperature, "K")


	plt.clf()
	plt.minorticks_on()
	plt.plot(time[time>T], voltage[time>T], ".")
	plt.xlabel("Zeit / s")

	plt.ylabel("gemessene Spannung / V")
	twin = secondary_y_axis(lambda x: x * 100)
	twin.set_ylabel("Signalspannung / uV")

	fit = Fit(LINEAR)
	fit.set_data(time, voltage, yerrors=1)
	fit = fit.filtered(time > T)
	fit.iterative_fit(1)
	fit.plot(plot_data=False)
	plt.xlim(T,time.max())
	plt.savefig("out/voltage_noise." + SAVETYPE)

	voltage = voltage - fit.eval(time)
	error_voltage = voltage[time > T].std()
	print("Voltage noise:", error_voltage, "V", "=", error_voltage*100, "uV")

	return error_temperature, error_voltage

def secondary_y_axis(func):
	axes = plt.gca()
	l, h = axes.get_ylim()
	twin = axes.twinx()
	tl, th = func(l), func(h)
	twin.set_ylim(tl, th)
	plt.sca(axes)
	return twin

def plot_cold_supra_imag(fit):
	# Imaginärteil Chi''
	temperature, voltage = load_samples("data/aufwaermen_supra_imag.dat")

	time = np.arange(0, len(temperature)) / 3

	plt.clf()
	plt.minorticks_on()
	plt.plot(temperature, voltage, '.')
	plt.xlabel("Temperatur / K")
	plt.ylabel("gemessene Spannung / V")
	twin = secondary_y_axis(lambda x: x * 100)
	twin.set_ylabel("Signalspannung / uV")
	plt.savefig("out/cold_supra_imag." + SAVETYPE)

	voltage -= fit.eval(temperature)

	F = np.logical_and(np.logical_and(temperature > 90, temperature < 100), voltage > 0.1)

	peakpos = ufloat((temperature[F]*voltage[F]).sum() / voltage[F].sum(), 0.2)
	print("Peak:", peakpos)

	plt.clf()
	plt.minorticks_on()
	plt.plot(temperature, voltage, '.', color="gray")
	plt.plot(temperature[F], voltage[F], '.', color="red")
	plt.axvline(peakpos.n, color="black", linestyle="--")
	plt.xlabel("Temperatur / K")
	plt.ylabel("gemessene Spannung / V")
	plt.xlim(80,120)
	twin = secondary_y_axis(lambda x: x * 100)
	twin.set_ylabel("Signalspannung / uV")
	info_box("Curietemperatur: ({:P}) K".format(peakpos), location="br")
	plt.savefig("out/cold_supra_imag_corrected." + SAVETYPE)

def plot_cold_supra_real(error_temperature, error_voltage):
	# Imaginärteil Chi''
	temperature, voltage = load_samples("data/aufwaermen_supra_real.dat")

	time = np.arange(0, len(temperature)) / 3
	plt.clf()
	plt.minorticks_on()
	plt.plot(temperature, voltage, '.')
	plt.xlabel("Temperatur / K")
	plt.ylabel("gemessene Spannung / V")
	twin = secondary_y_axis(lambda x: x * 100)
	twin.set_ylabel("Signalspannung / uV")
	plt.savefig("out/cold_supra_real." + SAVETYPE)

	plt.clf()

	F = time>10
	temperature = temperature[F]
	voltage = voltage[F]
	time = time[F]

	plt.clf()
	plt.minorticks_on()

	low = voltage[:10].mean()
	error_low = voltage[:10].std()

	chi = - voltage / low
	error_chi = np.sqrt(np.power(error_voltage / low, 2) + np.power(voltage / np.power(low, 2) * error_low, 2))

	plt.clf()
	plt.minorticks_on()
	plt.plot(temperature, chi, '.')
	plt.xlabel("Temperatur / K")
	plt.ylabel("Chi")
	plt.axhline(-1, color="black", linestyle="--")
	plt.axhline(-.5, color="black", linestyle="--")
	plt.axhline(0, color="black", linestyle="--")
	twin = secondary_y_axis(lambda x: -x*low*100)
	twin.set_ylabel("Signalspannung / uV")
	info_box("Chi = -1 bei {:.2f} uV".format(low*100), location="br")
	plt.savefig("out/cold_supra_real_chi." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	fit = Fit(LINEAR)
	fit.set_data(temperature, chi, xerrors=error_temperature, yerrors=error_chi)
	fit.set_labels(xlabel="Temperatur / K", ylabel="Chi")
	fit = fit.filtered(np.logical_and(temperature > 93, temperature < 94.5))
	fit.iterative_fit(5)
	plt.errorbar(temperature, chi, xerr=error_temperature, yerr=error_chi, fmt=',', color="gray")
	fit.plot(range=(92, 95.5), box="tl", plot_fit=True, plot_data=True, color="red")
	plt.xlim(90, 97.5)
	plt.axhline(-1, color="black", linestyle="--")
	plt.axhline(-.5, color="black", linestyle="--")
	plt.axhline(0, color="black", linestyle="--")
	T_c = (-0.5 - fit.uvalue("offset")) / fit.uvalue("slope")
	info_box("Curietemperatur: ({:P}) K".format(T_c), location="br")
	plt.savefig("out/cold_supra_real_chi_fit." + SAVETYPE)
	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tl")
	plt.savefig("out/cold_supra_real_chi_residual." + SAVETYPE)

	print("Curie Temperature:", T_c, "K")

	return ufloat(low, error_low)

def plot_sample(low, error_temperature, error_voltage):
	temperature, voltage = load_samples("data/aufwaermen_probe.dat")
	time = np.arange(0, len(temperature)) / 3

	plt.clf()
	plt.minorticks_on()
	plt.plot(temperature, voltage, '.')
	plt.xlabel("Temperatur / K")
	plt.ylabel("gemessene Spannung / V")
	twin = secondary_y_axis(lambda x: x * 100)
	twin.set_ylabel("Signalspannung / uV")
	plt.savefig("out/cold_sample_real." + SAVETYPE)

	chi = - voltage / low.n
	error_chi = np.sqrt(np.power(error_voltage / low.n, 2) + np.power(voltage / np.power(low.n, 2) * low.s, 2))

	plt.clf()
	plt.minorticks_on()
	plt.plot(temperature, chi, '.')
	plt.xlabel("Temperatur / K")
	plt.ylabel("Chi")
	#plt.axhline(0, color="black", linestyle="--")
	twin = secondary_y_axis(lambda x: -x*low.n*100)
	twin.set_ylabel("Signalspannung / uV")
	plt.savefig("out/cold_sample_real_chi." + SAVETYPE)

	plt.clf()
	plt.minorticks_on()
	plt.plot(temperature, 1/chi, '.')
	plt.xlabel("Temperatur / K")
	plt.ylabel("1 / Chi")
	plt.savefig("out/cold_sample_real_chiinv." + SAVETYPE)

	fit = Fit(LINEAR)
	fit.set_data(temperature, 1/chi, xerrors=error_temperature, yerrors=(np.power(1/chi, 2)*error_chi))
	fit.set_labels(xlabel="Temperatur / K", ylabel="1 / Chi")
	fit = fit.filtered(np.logical_and(temperature > 153, temperature < 165))
	fit.iterative_fit(5)
	fit.plot(range=(145, 180), box="tl", plot_data=True, plot_fit=True, color="red")

	zero = ufloat((1/chi[temperature < 140]).mean(), (1/chi[temperature < 140]).std())
	T_c = (zero-fit.uvalue("offset"))/fit.uvalue("slope")
	plt.axhline(zero.n, color="black", linestyle="--")
	plt.axvline(T_c.n, color="black", linestyle="--")

	plt.savefig("out/cold_sample_real_chiinv_fit." + SAVETYPE)
	plt.clf()
	plt.minorticks_on()
	fit.plot_residual(box="tl")
	plt.savefig("out/cold_sample_real_chiinv_residual." + SAVETYPE)

	T_c = ufloat(140,0.1)
	print("Curie temperature:", T_c, "K")

	temperature_difference = temperature-T_c.n
	error_temperature_difference = np.sqrt(np.power(error_temperature, 2) + np.power(T_c.s, 2))

	plt.clf()
	plt.minorticks_on()
	plt.plot(np.log(temperature_difference), np.log(chi), color="gray")
	plt.xlabel("ln(T-T_c)")
	plt.ylabel("ln(Chi)")
	plt.show()
	return

	fit = Fit(LINEAR)
	fit.set_data(np.log(temperature_difference), np.log(chi), xerrors=np.abs(error_temperature_difference / temperature_difference), yerrors=np.abs(error_chi / chi))
	fit.set_labels(xlabel="ln(T-T_c)", ylabel="ln(Chi)")
	fit = fit.filtered(np.logical_and(np.log(temperature_difference) > 1.5, np.log(temperature_difference) < 3))
	fit.iterative_fit(5)
	fit.plot(box="bl", color="red")
	plt.show()

if __name__=="__main__":
	#plot_teil1()
	#plot_tiefpass()
	#plot_hochpass()
	#plot_bandpass()
	plot_integration_time()
	#faktor()

	#error_temperature, error_voltage  = plot_cold_noise()
	#fit = plot_cold_empty(error_temperature, error_voltage)
	#plot_cold_supra_imag(fit)
	#low = plot_cold_supra_real(error_temperature, error_voltage)

	#plot_sample(low, error_temperature, error_voltage)
