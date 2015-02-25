import numpy as np

def dft(signal, frequency=1):
	N = len(signal)
	n = np.arange(N)
	fft = (signal[n]*np.exp(-2*np.pi*1j/N * np.outer(n,n))).sum(1)
	f = n/N * frequency
	return f, fft

def idft(fft, frequency=1):
	N = len(fft)
	n = np.arange(N)
	signal = (fft[n]*np.exp(+2*np.pi*1j/N * np.outer(n,n))).sum(1)/N
	t = n / frequency
	return t, signal

def smooth(data, minfreq, samplefreq=1):
	f, fy = dft(data, samplefreq)
	fy[f>minfreq] = 0.0
	_, y = idft(fy, samplefreq)
	return y.real
