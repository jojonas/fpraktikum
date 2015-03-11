from decimal import Decimal, Context, InvalidOperation
import re

def exponent2latex(number):
	return re.sub(r'([0-9\-.+]+)[Ee]([0-9\-+]+)', r'\1{\\times}10^{\2}', number)


def formatNumber(number, *, precision=2):
	c = Context(prec=precision)
	n = c.create_decimal_from_float(number)
	return n.to_eng_string(context=c)

def formatNumberPair(value, *errors, precision=2, separator=" +- "):
	c = Context(prec=precision)

	largest_error = None
	for error in errors:
		e = c.create_decimal_from_float(error)
		if e.is_finite():
			if largest_error:
				largest_error = largest_error.max_mag(e)
			else:
				largest_error = e

	results = []

	value = Decimal(value)

	if largest_error and value.is_finite():
		value = value.quantize(largest_error)
	results.append(value)

	for error in errors:
		error = Decimal(error)
		if largest_error and error.is_finite():
			error = error.quantize(largest_error)
		results.append(error)
	return tuple((str(r) for r in results))

def formatQuantityLatex(value, stat, sys=None, unit="", parenthesis=True, math=True, suffix=True):
	if not sys:
		value, stat = formatNumberPair(value, stat)
		retval = value + r' \pm ' + stat + (r'_\textrm{stat}' if suffix else '')

	else:
		value, stat, sys = formatNumberPair(value, stat, sys)
		retval = value + r' \pm ' + stat + (r'_\textrm{stat}' if suffix else '') + r' \pm ' + sys + (r'_\textrm{sys}' if suffix else '')

	retval = exponent2latex(retval)

	if parenthesis:
		retval = '(' + retval + ')'

	if unit:
		retval += " " + r'\enskip \mathrm{' + unit + '}'

	if unit and not parenthesis:
		print("Warning: parenthesis turned off for a formatted number with unit.\n", file=sys.stderr)

	if math:
		retval = r'$' + retval + r'$'

	return retval

def formatQuantity(value, stat, sys=None, unit="", parenthesis=True):
	if not sys:
		value, stat = formatNumberPair(value, stat)
		retval = value + ' \u00B1 ' + stat

	else:
		value, stat, sys = formatNumberPair(value, stat, sys)
		retval = value + ' \u00B1 ' + stat + ' \u00B1 ' + sys

	if parenthesis:
		retval = '(' + retval + ')'

	if unit:
		retval += " " + unit

	if unit and not parenthesis:
		print("Warning: parenthesis turned off for a formatted number with unit.\n", file=sys.stderr)

	return retval
