import os

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import download_file
from astropy.modeling import models
from astropy.modeling.tabular import Tabular1D

__all__ = ['FilterGenerator']

data_path = os.path.join(os.path.dirname(__file__), 'data', 'fft.fits')

include_facilities = ['2MASS', 'SLOAN', 'Kepler', 'TESS', 'HST', 'JWST',
                      'LSST', 'Keck', 'WISE', 'WFIRST', 'Spitzer', 'GAIA']

include_photsys = ['Bessel', 'Johnson', 'Cousins']


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

        return Filter(wavelength * u.Angstrom, transmittance)

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

        model = fft_model()
        return Filter(wavelength * u.Angstrom, model(wavelength), model=model)

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
        return Filter(true_transmittance['Wavelength'].data.data * u.Angstrom,
                      true_transmittance['Transmission'].data.data)

    def download_all_links(self, cache=True):
        from bs4 import BeautifulSoup
        bs = BeautifulSoup(open(download_file('http://svo2.cab.inta-csic.es/'
                                              'theory/fps3/fps.php')).read())
        all_options = [b.attrs['value'] for b in bs.find_all('option')]
        instruments = all_options[1:all_options.index('', 1)]
        facilities = all_options[all_options.index('', 1) + 1:
                                 all_options.index('', 150)]
        photometric_system = all_options[all_options.index('', 150) + 1:]

        filters = []

        for facility in include_facilities:
            url = ('http://svo2.cab.inta-csic.es/theory/fps3/'
                   'fps.php?Facility={0}'
                   .format(facility.replace(" ", "%20")))
            table = Table.read(download_file(url, timeout=300,
                                             cache=cache),
                               format='votable')
            filters.append(table)

        for photsys in include_photsys:
            url = ('http://svo2.cab.inta-csic.es/theory/fps3/'
                   'fps.php?PhotSystem={0}'
                   .format(photsys.replace(" ", "%20")))
            table = Table.read(download_file(url, timeout=300,
                                             cache=cache), format='votable')
            filters.append(table)

        links = [[i for i in filt['TrasmissionCurve']]
                 for filt in filters]

        return links

    def download_all_tables(self, links, cache=True):
        tables = dict()
        for facility_links in links:
            for link in facility_links:
                name = link.decode().split('=')[1]
                if (name not in tables.keys() and
                    name not in ['Scorpio/Comet.CO+']):
                    path = download_file(link.decode().replace('+', '&#43;'),
                                         cache=cache)
                    tables[name] = Table.read(path)

        return tables

    def fft_table(self, tables):
        rows = dict()
        n_terms = 10

        for k, v in tables.items():
            wl, tr = v['Wavelength'], v['Transmission']

            #     ax[0].plot(wl, tr, label=k, color='C0')

            diff_wl = np.diff(wl)

            delta_lambda = np.nanmedian(diff_wl[diff_wl != 0])
            lambda_0 = wl.min()
            n_lambda = len(wl)

            # Create a simplified wavelength grid:
            simplified_wavelength = np.arange(lambda_0, (n_lambda + 1) *
                                              delta_lambda + lambda_0,
                                              delta_lambda)
            tr_max = tr.max()

            # Interpolate transmittance onto simplified wavelength grid:
            tr_interp = np.interp(simplified_wavelength, wl, tr)

            # Take the DFT of the interpolated transmittance curve
            fft = np.fft.fft(tr_interp)[:n_terms]

            # Save results in a dictionary
            row = [n_lambda, lambda_0, delta_lambda, tr_max] + fft.tolist()
            rows[k] = row

        filtered_table = Table(rows=[[r] + rows[r] for i, r in enumerate(rows)
                                     if len(rows[r]) == n_terms + 4])

        bt = fits.BinTableHDU(data=filtered_table)

        return bt


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
        model : `~astropy.modeling.models.Model`
            Astropy model for the transmittance curve
        """
        self.wavelength = wavelength
        self.transmittance = transmittance
        self.model = model

    @property
    def table(self):
        return Tabular1D(points=self.wavelength,
                         lookup_table=self.transmittance)
