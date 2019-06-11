import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import os

__all__ = ['FourierFilter']

data_path = os.path.join(os.path.dirname(__file__), 'data', 'fft.fits')

class FourierFilter(object):
    def __init__(self, path=data_path): 
        self.path = path
        self.bintable = fits.getdata(path)
        
    def available_filters(self):
        return self.bintable.names
        
    def reconstruct(self, identifier):
        filt = self.bintable[identifier]
        n_lambda, lambda_0, delta_lambda, tr_max = filt[:4].real
        fft = filt[4:]

        wavelength = np.arange(lambda_0, (n_lambda + 1) * delta_lambda + lambda_0, delta_lambda)

        ifft = np.fft.ifft(fft, n=len(wavelength))

        flux = (ifft.real - ifft.real.min()) * tr_max / ifft.real.ptp()

        return wavelength, flux
