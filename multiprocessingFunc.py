import numpy as np
from scipy.interpolate import interp1d

## For some inane reason, this only worked when it was outside the class.
## This is the function that is passed to all of the processes.
def threeD_parallelCalculation(outputIndex, file, inputDict=[],):
	## Load spectra
	spectra_Bscan = np.loadtxt(file)
	
	## Background Subtraction
	if (inputDict['backgroundSubtraction']['enabled']):
		spectra_Bscan = spectra_Bscan - inputDict['backgroundSubtraction']['file']
	
	## Linearisation
	if (inputDict['linearisation']['enabled']):
		lin = inputDict['linearisation']
		spectra_Bscan = spectra_Bscan[::, lin['cuts'][0] : lin['cuts'][0]+lin['cuts'][1]]
		f = interp1d(lin['files']['nonlin'], spectra_Bscan, fill_value = 'extrapolate')
		spectra_Bscan = f(lin['files']['lin'])
	
	## Dispersion Compensation 
	if (inputDict['dispersionCompensation']['enabled']):
		dispComp = inputDict['dispersionCompensation']
		spectra_Bscan = spectra_Bscan[::, dispComp['cuts'][0] : dispComp['cuts'][0]+dispComp['cuts'][1]] * np.exp(-1j*dispComp['file'])
	
	## Calculate Fourier Transform
	fft_data_one = np.fft.fft(spectra_Bscan, n=inputDict['fftres'])
	abs_FFT = np.absolute(fft_data_one)
	
	if (inputDict['onlyOne']):
		return abs_FFT[:, inputDict['pix']]
	else:
		return abs_FFT