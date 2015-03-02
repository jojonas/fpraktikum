from decimal import Decimal, Context, InvalidOperation

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

if __name__=="__main__":
	for number in (5, 2000, -0.0000121233, 193890.2, -123123123123123, 0.0002312):
		print(number, "=>", formatNumber(number))
	print(formatNumberPair(12345678.0, 0.23))
