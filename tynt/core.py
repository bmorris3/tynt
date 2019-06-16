import os

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import download_file
from astropy.modeling import models
from astropy.modeling.tabular import Tabular1D

__all__ = ['FilterGenerator', 'Filter']

data_path = os.path.join(os.path.dirname(__file__), 'data', 'fft.fits')


class FilterGenerator(object):
    """
    Astronomical filter object generator.
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

    def reconstruct(self, identifier, model=False):
        """
        Reconstruct an approximate filter transmittance curve for
        a given filter.

        Parameters
        ----------
        identifier : str
            Name of the filter. To see available filters, run
            `~tynt.Filter.available_filters()`
        model : bool
            Construct a composite astropy model which approximates the
            transmittance curve.

        Returns
        -------
        filt : `~tynt.Filter`
            Astronomical filter object.
        """
        filt = list(self.table.loc[identifier])[1:]
        n_lambda, lambda_0, delta_lambda, tr_max = filt[:4]
        fft = filt[4:]
        astropy_model = None

        wavelength = np.arange(lambda_0, (n_lambda + 1) *
                               delta_lambda + lambda_0,
                               delta_lambda)
        N = len(wavelength)

        ifft = np.fft.ifft(fft, n=len(wavelength))

        transmittance = ((ifft.real - ifft.real.min()) * tr_max /
                         ifft.real.ptp())

        if model:
            m = (np.sum([models.Sine1D(amplitude=fft[i].real / N,
                                       frequency=i / N, phase=1 / 4)
                         for i in range(len(fft))]) -
                 np.sum([models.Sine1D(amplitude=fft[i].imag / N,
                                       frequency=i / N)
                         for i in range(len(fft))]))

            @models.custom_model
            def fft_model(x):
                """
                Approximate Fourier reconstruction of an astronomical filter

                Parameters
                ----------
                x : `~np.ndarray`
                    Wavelength in Angstroms.

                Returns
                -------
                transmittance : `~np.ndarray`
                    Transmittance curve
                """
                mo = m((x - wavelength.min()) /
                       (wavelength[1] - wavelength[0]))
                return (mo - mo.min()) * tr_max / mo.ptp()

            astropy_model = fft_model()

        return Filter(wavelength * u.Angstrom, transmittance,
                      model=astropy_model)

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
        filt : `~tynt.Filter`
            Astronomical filter object.
        """
        path = download_file('http://svo2.cab.inta-csic.es/'
                             'theory/fps3/fps.php?ID={0}'.format(identifier))

        true_transmittance = Table.read(path, format='votable')
        return Filter(true_transmittance['Wavelength'].data.data * u.Angstrom,
                      true_transmittance['Transmission'].data.data)


class Filter(object):
    """
    Astronomical filter object.
    """
    def __init__(self, wavelength, transmittance, model=None):
        """
        Parameters
        ----------
        wavelength : `~numpy.ndarray`
            Wavelength array
        transmittance : `~numpy.ndarray`
            Transmittance array
        model : `~astropy.modeling.Model`
            Astropy model for the transmittance curve
        """
        self.wavelength = wavelength
        self.transmittance = transmittance
        self.model = model

    @property
    def table(self):
        return Tabular1D(points=self.wavelength,
                         lookup_table=self.transmittance)
