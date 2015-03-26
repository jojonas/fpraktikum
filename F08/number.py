from decimal import Decimal, Context, InvalidOperation
from uncertainties import ufloat

def formatNumber(number, *, precision=2):
	c = Context(prec=precision)
	n = c.create_decimal_from_float(number)
	return n.to_eng_string(context=c)

def formatQuantityLatex(value, stat, sys=None, unit="", math=True):
	if sys:
		raise NotImplementedError
	uvalue = ufloat(value, stat)
	return formatUFloatLatex(uvalue, unit, math)

def formatQuantity(value, stat, sys=None, unit=""):
	if sys:
		raise NotImplementedError
	uvalue = ufloat(value, stat)
	return formatUFloat(uvalue, unit)

def formatUFloat(uvalue, unit=""):
	retval = '{:P}'.format(uvalue)
	if unit:
		retval = "(" + retval + ") " + unit
	return retval

def formatUFloatLatex(uvalue, unit="", math=True):
	retval = '{:L}'.format(uvalue)
	if unit:
		retval += r'\enskip\mathrm{' + unit + '}'
	if math:
		retval = r'$' + retval + r'$'
	return retval
