import sys
from warnings import warn

import math
import functools
import inspect
from itertools import zip_longest
from collections import OrderedDict

import matplotlib.pyplot as plt
from scipy.special import chdtrc
import numpy as np

import number
from funcs import *

_epsilon = math.sqrt(np.finfo(float).eps)

def derivative(func, x, dx):
	return (func(x+dx/2) - func(x-dx/2)) / dx

def chi2(value, y=0, err=1):
	return np.abs((value - y)/err)

def total_chi2(value, y=0, err=1):
	return np.power((value - y)/err, 2).sum()

def _minuit_costs(func, x, y, err):
	def f1(a):
		return total_chi2(func(x, a), y, err)
	def f2(a,b):
		return total_chi2(func(x, a, b), y, err)
	def f3(a,b,c):
		return total_chi2(func(x, a, b, c), y, err)
	def f4(a,b,c,d):
		return total_chi2(func(x, a, b, c, d), y, err)
	def f5(a,b,c,d,e):
		return total_chi2(func(x, a, b, c, d, e), y, err)
	def f6(a,b,c,d,e,f):
		return total_chi2(func(x, a, b, c, d, e, f), y, err)
	def f7(a,b,c,d,e,f,g):
		return total_chi2(func(x, a, b, c, d, e, f, g), y, err)

	func_map = {1: f1, 2: f2, 3: f3, 4: f4, 5: f5, 6: f6, 7: f7}
	names_map = ["a", "b", "c", "d", "e", "f", "g"]
	n_args = len(inspect.signature(func).parameters) - 1

	try:
		f = func_map[n_args]
		n = names_map[:n_args]
	except KeyError:
		raise ValueError("Invalid number of arguments, must be between 1 and 7.")

	return n, f

def _minimize_costs(func, x, y, err):
	def f(args):
		return total_chi2(func(x, *args), y, err)
	return f

def _lmfit_costs(func, x, y, err):
	def f(parameters):
		args = tuple(p.value for p in parameters.values())
		return chi2(func(x, *args), y, err)
	return f

def info_box(text, location='tl', margin=15, **custom):
	options = {
		"bbox": dict(boxstyle='square', facecolor='white'),
		"fontsize": 16,
		"linespacing": 1.4,
	}
	options.update(custom)

	if location == 'tr':
		plt.gca().annotate(text, xy=(1, 1), xytext=(-margin, -margin),
			xycoords='axes fraction', textcoords='offset points',
			horizontalalignment='right', verticalalignment='top', **options)
	elif location == 'br':
		plt.gca().annotate(text, xy=(1, 0), xytext=(-margin, margin),
			xycoords='axes fraction', textcoords='offset points',
			horizontalalignment='right', verticalalignment='bottom', **options)
	elif location == 'bl':
		plt.gca().annotate(text, xy=(0, 0), xytext=(margin, margin),
			xycoords='axes fraction', textcoords='offset points',
			horizontalalignment='left', verticalalignment='bottom', **options)
	else:
		plt.gca().annotate(text, xy=(0, 1), xytext=(margin, -margin),
			xycoords='axes fraction', textcoords='offset points',
			horizontalalignment='left', verticalalignment='top', **options)

class Fit:
	ERROR_PREFIX = "error_"

	def __init__(self, func):
		self._func = func
		self._params = OrderedDict()
		self._errors = OrderedDict()
		for i, (name, param) in enumerate(inspect.signature(func).parameters.items()):
			if i > 0:
				if param.default != inspect.Parameter.empty:
					self._params[name] = param.default
				else:
					self._params[name] = 0
				self._errors[self.ERROR_PREFIX + name] = 0
		self._chi2 = 0
		self._ndf = 0
		self._data_min = None
		self._data_max = None

	def __call__(self, x):
		return self.apply(x)

	def __setattr__(self, name, value):
		if name.startswith("_"):
			self.__dict__[name] = value
		elif name.startswith(self.ERROR_PREFIX) and name in self._errors:
			self._errors[name] = value
		elif name in self._params:
			self._params[name] = value

	def __getattr__(self, name):
		if name.startswith("_"):
			return self.__dict__[name]
		elif name.startswith(self.ERROR_PREFIX) and name in self._errors:
			return self._errors[name]
		elif name in self._params:
			return self._params[name]

	def apply(self, x):
		return self._func(x, *self.param_values())

	def param_values(self):
		return tuple(self._params.values())

	def param_names(self):
		return tuple(self._params.keys())

	def calculate_chi2(self, xdata, ydata, errors):
		return np.power((self.apply(xdata) - ydata)/errors, 2).sum()

	@property
	def chi2(self):
		return self._chi2

	@property
	def pvalue(self):
		return chdtrc(self._ndf, self._chi2)

	@property
	def ndf(self):
		return self._ndf

	@property
	def chi2ndf(self):
		return self._chi2 / self._ndf

	def _check_solver_error(self, error):
		if error is None:
			warn("Solver returned no error.")
		elif not isinstance(error, (int, float, complex)):
			warn("Solver returned invalid error type.")
		elif math.isinf(error):
			warn("Solver returned infinite error.")


	def fit(self, xdata, ydata, sigma, *, solver="minuit", method="migrad", **kwargs):
		assert xdata.shape == ydata.shape, "X and Y must have the same dimensions."
		assert sigma is None or sigma.shape == ydata.shape, "Errors and Y data must have the same dimensions."
		assert len(xdata) >= len(self._params), "Must have more data than free parameters."

		self._ndf = max(0, len(xdata) - len(self._params))

		if solver=="curve_fit":
			from scipy.optimize import curve_fit

			popt, pcov = curve_fit(self._func, xdata=xdata, ydata=ydata, p0=self.param_values(), sigma=sigma, absolute_sigma=True, **kwargs)
			for name, value, error in zip(self.param_names(), popt, np.sqrt(np.diag(pcov))):
				self._check_solver_error(error)

				self._params[name] = value
				self._errors[self.ERROR_PREFIX + name] = error

		elif solver=="minuit":
			import minuit

			minuit_names, costs = _minuit_costs(self._func, xdata, ydata, sigma)

			initial = {}
			for name, minuit_name in zip(self.param_names(), minuit_names):
				initial[minuit_name] = self._params[name]
				initial["error_" + minuit_name] = self._errors[self.ERROR_PREFIX + name]

			minuit = minuit.Minuit(costs, **initial)

			minuit.up = 1
			minuit.strategy = 1 # 0 = fast, 1 = default, 2 = thorough
			minuit.printMode = 0

			if method=="migrad":
				minuit.migrad(**kwargs)
				minuit.hesse()
			elif method=="simplex":
				minuit.simplex(**kwargs)
			elif method=="hesse":
				minuit.hesse(**kwargs)
			elif method=="minos":
				minuit.minos(**kwargs)
			else:
				raise ValueError("Invalid method.")

			for name, minuit_name in zip(self.param_names(), minuit_names):
				value = minuit.values[minuit_name]
				error = minuit.errors[minuit_name]

				self._check_solver_error(error)

				self._params[name] = value
				self._errors[self.ERROR_PREFIX + name] = error

		elif solver=="minimize":
			from scipy.optimize import minimize, OptimizeWarning

			costs = _minimize_costs(self._func, xdata, ydata, sigma)
			result = minimize(costs, x0=self.param_values(), method=method, jac=False, **kwargs)

			if not result.success:
				warn(result.message, OptimizeWarning)

			cov = np.linalg.inv(np.dot(result.jac.T, result.jac))

			for name, value, error in zip(self.param_names(), result.itervalues(), np.sqrt(np.diag(cov))):
				self._params[name] = value
				self._errors[self.ERROR_PREFIX + name] = error

		elif solver=="lmfit":
			from lmfit import minimize, fit_report, Parameters

			costs = _lmfit_costs(self._func, xdata, ydata, sigma)

			parameters = Parameters()
			for name, value in self._params.items():
				parameters.add(name, value=value, vary=True)

			result = minimize(costs, parameters, method=method)

			if not result.success:
				warn(result.message)

			for name in self.param_names():
				self._check_solver_error(parameters[name].stderr)

				self._params[name] = parameters[name].value
				self._errors[self.ERROR_PREFIX + name] = parameters[name].stderr


		else:
			raise ValueError("Invalid solver.")

		self._data_min = xdata.min()
		self._data_max = xdata.max()

		self._chi2 = self.calculate_chi2(xdata, ydata, sigma)

	def combine_errors(self, x, xerr, yerr, dx=_epsilon):
		slope = derivative(self.apply, x, dx)
		return np.sqrt(np.power(yerr,2) + np.power(xerr*slope, 2))

	def plot_residual(self, xdata, ydata, sigma=None, *, box=False, **plotopts):
		y = ydata - self.apply(xdata)
		plt.errorbar(xdata, y, yerr=sigma, **plotopts)
		plt.axhline(0, color="gray")
		if box:
			text =  r"$ \chi^2 $ / ndf = " + number.formatNumber(self.chi2) + " / " \
				+ number.formatNumber(self.ndf)
			text += "\n" + "p = " + number.formatNumber(self.pvalue)
			info_box(text, location=box)

	def plot(self, xmin=None, xmax=None, N=100, *, box=False, units={}, factors={}, **plotopts):
		if xmin is None:
			xmin = self._data_min
		if xmax is None:
			xmax = self._data_max

		x = np.linspace(xmin, xmax, N)
		y = self.apply(x)
		assert x.shape == y.shape
		plt.plot(x, y, '-', **plotopts)

		if box:
			lines = []
			for name, value in self._params.items():
				error = self._errors[self.ERROR_PREFIX + name]
				unit = units.get(name, "")
				factor = factors.get(name, 1)
				line = name + " = " + number.formatQuantityLatex(value*factor, error*factor, unit=unit)
				lines.append(line)
			text =  "\n".join(lines)
			info_box(text, location=box)

	def __str__(self):
		lines = []
		for name in self._params:
			value = self._params[name]
			error = self._errors[self.ERROR_PREFIX + name]
			lines.append(name + " = " + number.formatQuantity(value, error))
		return "<Fit '" + self._func.__name__ + "':: " + ", ".join(lines) + " >"


def linear_fit(x, y, yerr, xerr=None, N=10, slope=1, offset=0):
	fit = Fit(LINEAR)
	fit.slope = slope
	fit.offset = offset
	if xerr is None:
		fit.fit(x, y, yerr)
	else:
		for i in range(N):
			err = fit.combine_errors(x, xerr, yerr)
			if i < min(N/3, 3):
				fit.fit(x, y, err, solver="curve_fit")
			else:
				fit.fit(x, y, err, solver="minuit", method="migrad")

	return fit


def _simple_peak(data, index, sigma):
	lower = index-sigma
	upper = index+sigma

	maxindex = np.argmax(data[lower:upper])
	return maxindex + lower

def local_gauss_fit(data, mu=0, sigma=10, errors=None, N=5, f=1.0):
	fit = Fit(GAUSS)
	fit.mu = _simple_peak(data, mu, sigma)
	fit.sigma = sigma
	fit.A = data[fit.mu]

	for i in range(N):
		lower = max(int(fit.mu-f*fit.sigma), 0)
		upper = min(int(fit.mu+f*fit.sigma), len(data))

		xdata = np.arange(lower, upper)
		ydata = data[lower:upper]
		yerrors = errors[lower:upper]
		fit.fit(xdata, ydata, yerrors)
		fit.sigma = abs(fit.sigma)

	return fit
