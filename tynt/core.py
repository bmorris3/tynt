import os

import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import download_file
from astropy.modeling import models
from astropy.modeling.tabular import Tabular1D

__all__ = ['FilterGenerator', 'Filter']

data_path = os.path.join(os.path.dirname(__file__), 'data', 'fft.fits.zip')


class FilterGenerator:
    """
    Astronomical filter object generator.
    """
    def __init__(self, path=None):
        """
        Parameters
        ----------
        path : str
            Path to ``fft.fits`` file. If None, uses the package default data path.
        """
        if path is None:
            path = data_path

        self.path = path
        self.table = Table(fits.getdata(path))
        self.table.add_index('Filter name')

    def available_filters(self, as_dict=False, as_dataframe=False):
        """
        Return the available filters in the archive.

        Default returns a :py:class:`~numpy.ndarray`.

        Parameters
        ----------
        as_dict : bool
            Return available filters in a nested dictionary.
        as_dataframe : bool
            Return available filters in a pandas DataFrame.
        """
        if as_dataframe:
            as_dict = False

        if as_dict or as_dataframe:
            result = dict()
            for identifier in self.table['Filter name'].data:
                if '/' in identifier:
                    pre_slash, post_slash = identifier.split('/')
                    if pre_slash not in result:
                        result[pre_slash] = dict()
                    pre_period, post_period = (
                        post_slash.split('.')
                        if '.' in post_slash else (None, None)
                    )
                    if pre_period not in result[pre_slash]:
                        result[pre_slash].update({pre_period: []})
                    result[pre_slash][pre_period].append(post_period)
                else:
                    if 'OTHER' not in result:
                        result['OTHER'] = []
                    result['OTHER'].append(identifier)

            if as_dict:
                return result

            # otherwise, prep the dataframe:
            import pandas as pd
            reshaped = []
            for outerKey, innerDict in result.items():
                for innerKey, values in innerDict.items():
                    reshaped.append([outerKey, innerKey, ', '.join(values)])

            df = pd.DataFrame(reshaped, columns='group Set Filters'.split())
            df = df.set_index('Set')
            gr = df.groupby('group', group_keys=True)
            groups = dict()
            for name, group in gr:
                groups[name] = group.reindex(index=group.index)
            groups = pd.concat(groups)
            groups = groups.drop(columns=['group'])

            return groups

        return self.table['Filter name'].data

    def reconstruct(self, identifier, model=False):
        """
        Reconstruct an approximate filter transmittance curve for
        a given filter.

        Parameters
        ----------
        identifier : str
            Name of the filter. To see available filters, run
            :py:meth:`~tynt.FilterGenerator.available_filters`.
        model : bool
            Construct a composite astropy model which approximates the
            transmittance curve.

        Returns
        -------
        filt : ~tynt.Filter
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
                         np.ptp(ifft.real))

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
                x : ~np.ndarray
                    Wavelength in Angstroms.

                Returns
                -------
                transmittance : np.ndarray`
                    Transmittance curve
                """
                mo = m((x - wavelength.min()) /
                       (wavelength[1] - wavelength[0]))
                return (mo - mo.min()) * tr_max / np.ptp(mo)

            astropy_model = fft_model()

        return Filter(
            wavelength * u.Angstrom, transmittance, model=astropy_model
        )

    def download_true_transmittance(self, identifier, **kwargs):
        """
        Query the SVO service for a given filter,
        return the true transmittance curve as a
        :py:class:`~tynt.Filter`.

        Parameters
        ----------
        identifier : str
            Name of the filter. To see available filters, run
            :py:meth:`~tynt.FilterGenerator.available_filters`
        **kwargs : dict
            Passed to :py:func:`astropy.utils.data.download_file`
        """
        path = download_file('http://svo2.cab.inta-csic.es/'
                             f'theory/fps3/fps.php?ID={identifier}', **kwargs)

        true_transmittance = Table.read(path, format='votable')
        return Filter(
            true_transmittance['Wavelength'].data.data * u.Angstrom,
            true_transmittance['Transmission'].data.data
        )


class Filter:
    """
    Astronomical filter.
    """
    @u.quantity_input(wavelength=u.m)
    def __init__(self, wavelength, transmittance, model=None):
        """
        Parameters
        ----------
        wavelength : ~astropy.units.Quantity
            Wavelength array
        transmittance : ~numpy.ndarray
            Transmittance array
        model : ~astropy.modeling.Model
            Astropy model for the transmittance curve
        """
        self.wavelength = wavelength
        self.transmittance = transmittance
        self.model = model

    @property
    def table(self):
        """
        Filter transmittance represented in the astropy modeling framework.

        Returns
        -------
        tab : ~astropy.modeling.tabular.Tabular1D
            Lookup table for astropy modeling.
        """
        return Tabular1D(points=self.wavelength,
                         lookup_table=self.transmittance)
