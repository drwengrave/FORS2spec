import os
os.chdir('output_GRIS_300V/')
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
from matplotlib.colors import LogNorm
from pylab import figure, cm

from scipy import signal
from scipy.ndimage import median_filter
from numpy.polynomial import Chebyshev

#load 2D image
hdu = fits.open('mapped_flux_sci_1_GRIS_300V.fits')
hduerr = fits.open('mapped_flux_sci_1_GRIS_300V.fits')

#hdu = fits.open('VISOB2_binned.fits')
data = hdu[0].data
error = hduerr[1].data
head = hdu[0].header

valmean = np.abs(np.mean(data))
valmax = np.max(data)

plt.figure()
plt.imshow(data, cmap=cm.gray, norm=LogNorm(vmin=-0.01, vmax=1))
plt.show()

#trace of spectrum
trace = data[108:121,:]
error_trace = error[108:121,:]
bkgd = data[150:163]
bkgm = np.median(bkgd, axis=0)
vtrace = data[:,:]
vmean = np.mean(vtrace, axis=1)
bmean = np.mean(bkgd, axis=1)

#plot vtrace
plt.figure(figsize=(8,6))
plt.plot(np.arange(len(vmean)),vmean)
plt.show()

#flux
fluxU = np.sum(trace, axis=0)
errfluxU = np.sqrt(np.sum(error_trace**2, axis=0))
#fluxUb = np.sum(trace-bkgm, axis=0)
vfluxU = np.mean(vtrace, axis=1)

#wavelength
waveU = []
for i in range(len(fluxU)):
    waveU.append(hdu[0].header['CRVAL1'] + i*hdu[0].header['CD1_1'])

#plot
plt.figure(figsize=(8,6))
plt.plot(waveU,fluxU,lw=0.7,label='flux')
plt.plot(waveU,errfluxU,lw=0.7,label='error')
#plt.plot(np.arange(len(vfluxU)),vfluxU)
#plt.plot(waveU[1200:],fluxUb[1200:]+0.5e-16)
#plt.xlim(5318,8422)
plt.ylim(-0.01,1.79)
plt.ylabel(r'Flux')
plt.xlabel(r'Wavelength ($\AA$)')
plt.legend()
plt.show()

#export
f = open('spec_1D.txt', 'w')
for i in range(len(fluxU)):
    f.write('{:.5e}    {:.5e}    {:.5e}\n'.format(waveU[i], fluxU[i], errfluxU[i]))

f.close()
