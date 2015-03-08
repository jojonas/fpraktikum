import sys
import re

import math
import functools
import inspect
from itertools import zip_longest
from collections import OrderedDict

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import chdtrc
import numpy as np

import number
from funcs import *

def derivative(func, x, dx):
	return (func(x+dx/2) - func(x-dx/2)) / dx

def exponent2latex(number):
	return re.sub(r'([0-9\-.+]+)[Ee]([0-9\-+]+)', r'\1{\\times}10^{\2}', number)

def format_error(value, stat, sys=None, unit="", parenthesis=True, surroundmath=True):
	if unit:
		unit = r'\enskip \mathrm{' + unit + '}'
	if not sys:
		value, stat = number.formatNumberPair(value, stat)
		retval = value + r' \pm ' + stat + r'_\textrm{stat}'
	else:
		value, stat, sys = number.formatNumberPair(value, stat, sys)
		retval = value + r' \pm ' + stat + r'_\textrm{stat} \pm ' + sys + r'_\textrm{sys}'

	retval = exponent2latex(retval)

	if parenthesis:
		retval = '(' + retval + ')'

	if unit:
		retval += " " + unit

	if unit and not parenthesis:
		print("Warning: parenthesis turned off for a formatted number with unit.\n", file=sys.stderr)

	if surroundmath:
		retval = r'$' + retval + r'$'

	return retval

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
	ERROR_PREFIX = "sigma_"

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
		self._chisq = 0
		self._ndf = 0

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
		return self._func(x, *self.paramtuple())

	def paramtuple(self):
		return tuple(self._params.values())

	@property
	def chisq(self):
		return self._chisq

	@property
	def ndf(self):
		return self._ndf

	@property
	def chisqndf(self):
		return self._chisq / self._ndf

	def fit(self, xdata, ydata, sigma):
		assert xdata.shape == ydata.shape, "X and Y must have the same dimensions."
		assert sigma is None or sigma.shape == ydata.shape, "Errors and Y data must have the same dimensions."
		assert len(xdata) >= len(self._params), "Must have more data than free parameters."

		popt, pcov = curve_fit(self._func, xdata=xdata, ydata=ydata, p0=self.paramtuple(), sigma=sigma, absolute_sigma=True)
		for name, value, error in zip(self._params, popt, np.sqrt(np.diag(pcov))):
			self._params[name] = value
			if math.isinf(error):
				print("Warning: The error of '" + name + "' is infinite!\n" + \
				"\tA likely cause of this is that one of the fit parameters has no influence on the result.\n" + \
				"\tDouble-check your fitting function or try changing the inital parameters.\n"
					, file=sys.stderr)
			self._errors[self.ERROR_PREFIX + name] = error

		self._ndf = len(xdata) - len(self._params)
		self._chisq = np.power((self.apply(xdata) - ydata)/sigma, 2).sum()

	def combine_errors(self, x, xerr, yerr, dx=0.01):
		slope = derivative(self.apply, x, dx)
		return np.sqrt(np.power(yerr,2) + np.power(xerr*slope, 2))

	def plot_residuums(self, xdata, ydata, sigma=None, *, box=False, **plotopts):
		y = ydata - self.apply(xdata)
		plt.errorbar(xdata, y, yerr=sigma, **plotopts)
		plt.axhline(0, color="gray")
		if box:
			text =  r"$ \chi^2 $ / ndf = " + number.formatNumber(self.chisq) + " / " \
				+ number.formatNumber(self.ndf)
			text += "\n" + "p = " + number.formatNumber(chdtrc(self.ndf, self.chisq))
			info_box(text, location=box)

	def plot(self, xmin=-1, xmax=1, N=100, *, box=False, units={}, **plotopts):
		x = np.linspace(xmin, xmax, N)
		y = self.apply(x)
		plt.plot(x, y, '-', **plotopts)

		if box:
			lines = []
			for name, value in self._params.items():
				error = self._errors[self.ERROR_PREFIX + name]
				unit = units.get(name, "")
				line = name + " = " + format_error(value, error, unit=unit)
				lines.append(line)
			text =  "\n".join(lines)
			info_box(text, location=box)

def linear_fit(x, y, yerr, xerr=None, N=10):
	fit = Fit(LINEAR)
	if xerr is None:
		fit.fit(x, y, yerr)
	else:
		for i in range(N):
			err = fit.combine_errors(x, xerr, yerr)
			fit.fit(x, y, err)

	return fit


def _simple_peak(data, index, sigma):
	lower = index-sigma
	upper = index+sigma

	maxindex = np.argmax(data[lower:upper])
	return maxindex + lower

def local_gauss_fit(data, mu=0, sigma=10, A=100, xerrors=None, yerrors=None, N=5, f=1.0):
	fit = Fit(GAUSS)
	fit.mu = _simple_peak(data, mu, sigma)
	fit.sigma = sigma
	fit.A = A

	for i in range(N):
		lower = max(int(fit.mu-f*fit.sigma), 0)
		upper = min(int(fit.mu+f*fit.sigma), len(data))

		xdata = np.arange(lower, upper)
		ydata = data[lower:upper]
		errors = fit.combine_errors(xdata, xerrors[lower:upper], yerrors[lower:upper])
		fit.fit(xdata, ydata, errors)
		fit.sigma = abs(fit.sigma)

	return fit

if __name__=="__main__":
	print(exponent2latex("123e3, 1e-55"))
