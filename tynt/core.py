import os

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import download_file
from astropy.modeling import models


__all__ = ['Filter']

data_path = os.path.join(os.path.dirname(__file__), 'data', 'fft.fits')


class Filter(object):
    """
    Astronomical filter object.
    """

    def __init__(self, path=data_path):
        """
        Parameters
        ----------
        path : str (optional)
            Path to `fft.fits` file
        """
        self.path = path
        self.table = Table(fits.getdata(path))
        self.table.add_index('col0')

    def available_filters(self):
        """
        Return the available filters in the archive
        """
        return self.table['col0'].data

    def reconstruct(self, identifier):
        """
        Reconstruct an approximate filter transmittance curve for
        a given filter.

        Parameters
        ----------
        identifier : str
            Name of the filter. To see available filters, run
            `~tynt.Filter.available_filters()`

        Returns
        -------
        wavelength : `~numpy.ndarray`
            Wavelength array in Angstroms
        transmittance : `~numpy.ndarray`
            Approximate transmittance as a function of wavelength
        """
        filt = list(self.table.loc[identifier])[1:]
        n_lambda, lambda_0, delta_lambda, tr_max = filt[:4]
        fft = filt[4:]

        wavelength = np.arange(lambda_0, (n_lambda + 1) *
                               delta_lambda + lambda_0,
                               delta_lambda)

        ifft = np.fft.ifft(fft, n=len(wavelength))

        transmittance = ((ifft.real - ifft.real.min()) * tr_max /
                         ifft.real.ptp())

        return wavelength, transmittance

    def model(self, identifier):
        """
        Reconstruct an approximate filter transmittance curve using
        astropy models for a given filter.

        Parameters
        ----------
        identifier : str
            Name of the filter. To see available filters, run
            `~tynt.Filter.available_filters()`

        Returns
        -------
        wavelength : `~numpy.ndarray`
            Wavelength array in Angstroms
        model : `~astropy.modeling.Model`
            Approximate astropy model representing the
            filter transmission curve
        """
        filt = list(self.table.loc[identifier])[1:]
        n_lambda, lambda_0, delta_lambda, tr_max = filt[:4]
        fft = filt[4:]

        wavelength = np.arange(lambda_0, (n_lambda + 1) *
                               delta_lambda + lambda_0,
                               delta_lambda)
        N = len(wavelength)

        m = (np.sum([models.Sine1D(amplitude=fft[i].real/N,
                                   frequency=i/N, phase=1/4)
                     for i in range(len(fft))]) -
             np.sum([models.Sine1D(amplitude=fft[i].imag/N,
                                   frequency=i/N)
                     for i in range(len(fft))]))

        @models.custom_model
        def fft_model(x):
            mo = m((x - wavelength.min()) /
                   (wavelength[1] - wavelength[0]))
            return (mo - mo.min()) * tr_max / mo.ptp()

        return wavelength, fft_model()

    def download_true_transmittance(self, identifier):
        """
        Query the SVO service for a given filter,
        return the true transmittance curve.

        Parameters
        ----------
        identifier : str
            Name of the filter. To see available filters, run
            `~tynt.Filter.available_filters()`

        Returns
        -------
        wavelength : `~numpy.ndarray`
            True wavelength array in Angstroms
        transmittance : `~numpy.ndarray`
            True transmittance as a function of wavelength
        """
        path = download_file('http://svo2.cab.inta-csic.es/'
                             'theory/fps3/fps.php?ID={0}'.format(identifier))

        true_transmittance = Table.read(path, format='votable')
        return (true_transmittance['Wavelength'].data.data,
                true_transmittance['Transmission'].data.data)
