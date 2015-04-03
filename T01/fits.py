import sys
from warnings import warn

import math
import copy
import functools
from itertools import zip_longest
from collections import OrderedDict
import inspect

import matplotlib.pyplot as plt
from scipy.special import chdtrc
import numpy as np

import number
from funcs import *

_epsilon = math.sqrt(np.finfo(float).eps)

def derivative(f, x, h):
	return (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h)) / (12*h)

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
	def __init__(self, func):
		assert callable(func), "Parameter 'func' must be a callable object."

		self._func = func
		self._param_values = OrderedDict()
		self._param_errors = OrderedDict()

		self._load_from_signature(func)

		self._x = None
		self._y = None
		self._xerrors = None
		self._yerrors = None

		self._xlabel = None
		self._ylabel = None

		self._fitted = False

	def __call__(self, x):
		return self.eval(x)

	def __len__(self):
		return len(self._x)

	def _load_from_signature(self, func):
		self._param_values.clear()
		self._param_errors.clear()

		for i, (name, param) in enumerate(inspect.signature(func).parameters.items()):
			if i > 0:
				if param.default != inspect.Parameter.empty:
					self._param_values[name] = param.default
				else:
					self._param_values[name] = 0
				self._param_errors[name] = 0

	def _check_solver_value(self, value):
		if value is None:
			warn("Solver returned no value.")
		elif not isinstance(value, (int, float, complex)):
			warn("Solver returned unknon value type.")
		elif math.isnan(value):
			warn("Solver returned 'not a number'-value.")
		elif math.isinf(value):
			warn("Solver returned infinite value.")

	def _assert_data_sanity(self):
		assert callable(self._func), "No function specified."
		assert len(self._param_values) == len(self._param_errors), "Parameter count {:d} and error count {:d} are not equal.".format(len(self._param_values), len(self._param_errors))
		assert self._xlabel is None or isinstance(self._xlabel, str), "X label must be a string."
		assert self._ylabel is None or isinstance(self._ylabel, str), "Y label must be a string."
		if self._x is not None or self._y is not None:
			assert self._x is not None and self._y is not None, "If X is set, Y must be set too (and vice-versa)."
			assert isinstance(self._x, np.ndarray), "X data must be a numpy array."
			assert isinstance(self._y, np.ndarray), "Y data must be a numpy array."
			assert self._x.shape == self._y.shape, "X and Y must have the same dimensions."
			assert len(self._x) >= len(self._param_values), "Must have more data than free parameters."
		if self._xerrors is not None:
			assert isinstance(self._xerrors, np.ndarray), "X errors must be a numpy array."
			assert self._xerrors.shape == self._x.shape, "Errors and X data must have the same dimensions."
		if self._yerrors is not None:
			assert isinstance(self._yerrors, np.ndarray), "Y errors must be a numpy array."
			assert self._yerrors.shape == self._yerrors.shape, "Errors and Y data must have the same dimensions."

	def __getitem__(self, value):
		return self.filtered(where=value)

	def filtered(self, where):
		new = self.__class__(self._func)
		new._param_values = copy.deepcopy(self._param_values)
		new._param_errors = copy.deepcopy(self._param_errors)
		new._xlabel = self._xlabel
		new._ylabel = self._ylabel

		# slicing only! (keeps references!)
		if self._x is not None:
			new._x = self._x[where]
		if self._y is not None:
			new._y = self._y[where]
		if self._xerrors is not None:
			new._xerrors = self._xerrors[where]
		if self._yerrors is not None:
			new._yerrors = self._yerrors[where]
		return new

	def eval(self, x):
		return self._func(x, *self.param_values())

	def ueval(self, x):
		return self._func(x, *self.param_uvalues())

	def diff(self, x, dx=_epsilon):
		return derivative(self.eval, x, dx)

	def udiff(self, x, dx=_epsilon):
		return derivative(self.ueval, x, dx)

	def param_values(self):
		return tuple(self._param_values.values())

	def param_errors(self):
		return tuple(self._param_errors.values())

	def param_names(self):
		return tuple(self._param_values.keys())

	def param_uvalues(self):
		from uncertainties import ufloat
		return tuple((ufloat(v,e) for v,e in zip(self.param_values(), self.param_errors())))

	def value(self, name):
		return self._param_values[name]

	def error(self, name):
		return self._param_errors[name]

	def uvalue(self, name):
		from uncertainties import ufloat
		return ufloat(self._param_values[name], self._param_errors[name])

	@property
	def chi2(self):
		assert self._x is not None and self._y is not None, "Cannot compute Chi^2 without data."
		if not self._fitted:
			warn("Calculating Chi^2 of non-fit.")

		errors = self.combined_errors(allow_zero=False)
		return np.power((self.eval(self._x) - self._y)/errors, 2).sum()

	@property
	def pvalue(self):
		return chdtrc(self.ndf, self.chi2)

	@property
	def ndf(self):
		if self._x is None:
			return 0
		else:
			return max(0, len(self._x) - len(self._param_values))

	@property
	def chi2ndf(self):
		return self.chi2 / self.ndf

	@property
	def usetex(self):
		return plt.rcParams["text.usetex"]


	def set_data(self, xdata=None, ydata=None, *, xerrors=None, yerrors=None):
		ARRAY_TYPE = (np.ndarray, list, tuple)
		NUMBER_TYPE = (float, int)

		if isinstance(xdata, ARRAY_TYPE):
			self._x = np.array(xdata)
		else:
			self._x = None

		if isinstance(ydata, ARRAY_TYPE):
			self._y = np.array(ydata)
		else:
			self._y = None

		if isinstance(xerrors, NUMBER_TYPE):
			self._xerrors = np.ones_like(xdata) * xerrors
		elif isinstance(xerrors, ARRAY_TYPE):
			self._xerrors = np.array(xerrors)
		else:
			self._xerrors = None

		if isinstance(yerrors, NUMBER_TYPE):
			self._yerrors = np.ones_like(ydata) * yerrors
		elif isinstance(yerrors, ARRAY_TYPE):
			self._yerrors = np.array(yerrors)
		else:
			self._yerrors = None


	def set_labels(self, xlabel='', ylabel=''):
		assert isinstance(xlabel, str) and isinstance(ylabel, str), "Labels must be strings!"

		self._xlabel = xlabel
		self._ylabel = ylabel

	def set_params(self, **kwargs):
		for name, value in kwargs.items():
			assert name in self._param_values, "Unknown parameter '{:s}'".format(name)
			self._param_values[name] = value

	def set_errors(self, **kwargs):
		for name, value in kwargs.items():
			self._param_errors[name] = value

	def iterative_fit(self, N, **kwargs):
		for _ in range(N):
			self.fit(**kwargs)

	def fit(self, *, solver="minuit", method="migrad", **kwargs):
		assert self._x is not None and self._y is not None, "No data to fit."

		self._assert_data_sanity()

		errors = self.combined_errors(allow_zero=False)

		if solver=="curve_fit":
			from scipy.optimize import curve_fit

			popt, pcov = curve_fit(self._func, xdata=self._x, ydata=self._y, p0=self.param_values(), sigma=errors, absolute_sigma=True, **kwargs)
			for name, value, error in zip(self.param_names(), popt, np.sqrt(np.diag(pcov))):
				self._check_solver_value(value)
				self._check_solver_value(error)

				self._param_values[name] = value
				self._param_errors[name] = error

		elif solver=="minuit":
			import minuit

			minuit_names, costs = _minuit_costs(self._func, self._x, self._y, errors)

			initial = {}
			for name, minuit_name in zip(self.param_names(), minuit_names):
				initial[minuit_name] = self._param_values[name]
				initial["error_" + minuit_name] = self._param_errors[name]

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

				self._check_solver_value(value)
				self._check_solver_value(error)

				self._param_values[name] = value
				self._param_errors[name] = error

		elif solver=="minimize":
			from scipy.optimize import minimize, OptimizeWarning

			costs = _minimize_costs(self._func, self._x, self._y, errors)
			result = minimize(costs, x0=self.param_values(), method=method, jac=False, **kwargs)

			if not result.success:
				warn(result.message, OptimizeWarning)

			cov = np.linalg.inv(np.dot(result.jac.T, result.jac))

			for name, value, error in zip(self.param_names(), result.itervalues(), np.sqrt(np.diag(cov))):
				self._param_values[name] = value
				self._param_errors[name] = error

		elif solver=="lmfit":
			from lmfit import minimize, fit_report, Parameters

			costs = _lmfit_costs(self._func, self._x, self._y, errors)

			parameters = Parameters()
			for name, value in self._param_values.items():
				parameters.add(name, value=value, vary=True)

			result = minimize(costs, parameters, method=method)

			if not result.success:
				warn(result.message)

			for name in self.param_names():
				self._check_solver_value(parameters[name].value)
				self._check_solver_value(parameters[name].stderr)

				self._param_values[name] = parameters[name].value
				self._param_errors[name] = parameters[name].stderr

		else:
			raise ValueError("Invalid solver.")

		self._fitted = True

	def combined_errors(self, allow_zero=True):
		assert self._x is not None, "Cannot combine errors without having X data."

		if self._xerrors is not None:
			xerr = self._xerrors
		else:
			xerr = 0
		if self._yerrors is not None:
			yerr = self._yerrors
		else:
			yerr = 0

		if not allow_zero and self._xerrors is None and self._yerrors is None:
			raise AssertionError("No errors present for combination.")

		slope = derivative(self.eval, self._x, _epsilon)
		return np.sqrt(np.power(yerr,2) + np.power(xerr*slope, 2))

	def plot_residual(self, *, box=False, **plotopts):
		assert self._y is not None and self._x is not None, "Cannot plot residual graph without data."

		if not self._fitted:
			warn("Plotting residuals of non-fit.")

		y = self._y - self.eval(self._x)
		errors = self.combined_errors(allow_zero=True)

		options = {
			"color": "black",
			"fmt": ".",
		}
		options.update(plotopts)

		plt.errorbar(self._x, y, yerr=errors, **options)
		plt.axhline(0, color="gray")
		if self._xlabel:
			plt.xlabel(self._xlabel)
		if self._ylabel:
			plt.ylabel(self._ylabel)
		if box:
			if self.usetex:
				text = r"$ \chi^2 $ / ndf = " + number.formatNumber(self.chi2) + " / " + number.formatNumber(self.ndf)
			else:
				text = "chi\u00B2 / ndf = " + number.formatNumber(self.chi2) + " / " + number.formatNumber(self.ndf)

			text += "\n" + "p = " + number.formatNumber(self.pvalue)
			info_box(text, location=box)

	def plot(self, range=None, N=100, *, plot_data=True, plot_fit=True, box=False, units={}, factors={}, **plotopts):
		assert plot_data or plot_fit, "Nothing to plot!"

		if plot_data and self._x is not None and self._y is not None:
			options = {
				"color": "black",
				"fmt": ","
			}
			options.update(plotopts)
			plt.errorbar(self._x, self._y, xerr=self._xerrors, yerr=self._yerrors, **options)

		if plot_fit:
			if range is None:
				assert self._x is not None, "Needs either data or range for plotting."
				xmin = self._x.min()
				xmax = self._x.max()
			else:
				xmin, xmax = range

			x = np.linspace(xmin, xmax, N)
			y = self.eval(x)
			assert x.shape == y.shape

			options = {
				"color": "red",
			}
			options.update(plotopts)
			plt.plot(x, y, '-', **options)

			if box:
				lines = []
				for name, value in self._param_values.items():
					error = self._param_errors[name]
					unit = units.get(name, "")
					factor = factors.get(name, 1)

					if self.usetex:
						line = name + " = " + number.formatQuantityLatex(value*factor, error*factor, unit=unit)
					else:
						line = name + " = " + number.formatQuantity(value*factor, error*factor, unit=unit)

					lines.append(line)
				text =  "\n".join(lines)
				info_box(text, location=box)

		if self._xlabel:
			plt.xlabel(self._xlabel)
		if self._ylabel:
			plt.ylabel(self._ylabel)

	def __str__(self):
		lines = []
		for name in self._param_values:
			lines.append(name + " = " + str(self.uvalue(name)))
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
