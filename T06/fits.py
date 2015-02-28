import math
import functools
import inspect
from itertools import zip_longest
from collections import OrderedDict

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

from funcs import *

def derivative(func, x, dx):
	return (func(x+dx/2) - func(x-dx/2)) / dx

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
		self._ndof = 0

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
	def ndof(self):
		return self._ndof

	@property
	def chisqndof(self):
		return self._chisq / self._ndof

	def fit(self, xdata, ydata, sigma=None):
		assert xdata.shape == ydata.shape, "X and Y must have the same dimensions."
		assert sigma is None or sigma.shape == ydata.shape, "Errors and Y data must have the same dimensions."
		assert len(xdata) >= len(self._params), "Must have more data than free parameters."

		popt, pcov = curve_fit(self._func, xdata=xdata, ydata=ydata, p0=self.paramtuple(), sigma=sigma, absolute_sigma=True)
		for name, value, error in zip(self._params, popt, np.sqrt(np.diag(pcov))):
			self._params[name] = value
			self._errors[self.ERROR_PREFIX + name] = error
		self._ndof = len(xdata) - len(self._params)
		self._chisq = np.power((self.apply(xdata) - ydata)/sigma, 2).sum()

	def combine_errors(self, x, xerr, yerr, dx=0.01):
		slope = derivative(self.apply, x, dx)
		return np.sqrt(np.power(yerr,2) + np.power(xerr*slope, 2))

	def plot_residuums(self, xdata, ydata, sigma=None, **plotopts):
		y = ydata - self.apply(xdata)
		plt.errorbar(xdata, y, yerr=sigma, **plotopts)
		plt.axhline(0, color="gray")

	def plot(self, xmin=-1, xmax=1, N=100, **plotopts):
		x = np.linspace(xmin, xmax, N)
		y = self.apply(x)
		plt.plot(x, y, '-', **plotopts)

		props = dict(boxstyle='square', facecolor='white')
		text = "\n".join((name + " = %.2f" % value for name, value in self._params.items()))
		plt.gca().text(0.05, 0.95, text, transform=plt.gca().transAxes, fontsize=14,
			verticalalignment='top', bbox=props)

def linear_fit(x, y, yerr, xerr=None, N=10):
	fit = Fit(LINEAR)
	if xerr is None:
		fit.fit(x, y, yerr)
	else:
		for i in range(N):
			err = fit.combine_errors(x, xerr, yerr)
			fit.fit(x, y, err)

	return fit


def _simple_peak(data, index):
	for i in range(len(data)):
		if index+1 < len(data) and data[index+1] > data[index]:
			index += 1
		elif index-1 >= 0 and data[index-1] > data[index]:
			index -= 1
		else:
			break
	return index

def local_gauss_fit(data, mu=0, sigma=10, A=100, N=5, errors=None, f=1.0):
	fit = Fit(GAUSS)
	fit.mu = _simple_peak(data, mu)
	fit.sigma = sigma
	fit.A = A

	for i in range(N):
		lower = max(int(fit.mu-f*fit.sigma), 0)
		upper = min(int(fit.mu+f*fit.sigma), len(data))

		xdata = np.arange(lower, upper)
		ydata = data[lower:upper]
		yerrors = errors[lower:upper]
		fit.fit(xdata, ydata, yerrors)
		fit.sigma = abs(fit.sigma)

	return fit
