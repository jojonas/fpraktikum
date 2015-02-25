import numpy as np

def linreg(x, y, err):
	x = np.array(x)
	y = np.array(y)
	err = np.array(err)
	var = np.power(err,2)
	
	#Sx = (x/var).sum()
	#S = (1/var).sum()
	#correction = Sx/S
	#x -= correction

	S = (1/var).sum()
	Sy = (y/var).sum()
	Sx = (x/var).sum()
	Sxy = (x*y/var).sum()
	Sxx = (x*x/var).sum()

	delta = S*Sxx-Sx*Sx

	slope = (S*Sxy-Sx*Sy) / delta
	offset = (Sxx*Sy-Sx*Sxy) / delta
	error_slope = np.sqrt(S/delta)
	error_offset = np.sqrt(Sxx/delta)

	#x += correction
	#offset = correction*slope

	return (slope, offset), (error_slope, error_offset)

def linreg2(x, y, xerr, yerr, N=10, slope=1, offset=0):
	x = np.array(x)
	y = np.array(y)
	xerr = np.array(xerr)
	yerr = np.array(yerr)

	for i in range(1):
		err = np.sqrt(np.power(xerr,2) + np.power(yerr*slope, 2))
		popt, perr = linreg(x, y, err)
		slope, offset = popt
	return popt, perr
