import os

import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.io import fits
from astropy.modeling import models
from astropy.modeling.tabular import Tabular1D
from astropy.table import Table
from astropy.utils.data import download_file
from astropy.visualization import quantity_support

__all__ = ['FilterGenerator', 'Filter']

data_path = os.path.join(os.path.dirname(__file__), 'data', 'fft.fits.zip')


class FilterGenerator:
    """
    Astronomical filter object generator.
    """
    def __init__(self, path=None, name=None):
        """
        Parameters
        ----------
        path : str
            Path to ``fft.fits`` file. If None, uses the package default data path.
        """
        if path is None:
            path = data_path

        self.path = path
        self.name = name
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

        filter_set, filter_name = identifier.split('/')
        return Filter(
            wavelength * u.Angstrom, transmittance,
            filter_set=filter_set, filter_name=identifier,
            model=astropy_model
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

        filter_set, filter_name = identifier.split('/')

        return Filter(
            true_transmittance['Wavelength'].data.data * u.Angstrom,
            true_transmittance['Transmission'].data.data,
            filter_set=filter_set, filter_name=identifier
        )


_filter_generator = FilterGenerator()


class Filter:
    """
    Astronomical filter.
    """
    @u.quantity_input(wavelength=u.m)
    def __init__(
            self, wavelength=None, transmittance=None,
            model=None, filter_set=None, filter_name=None
    ):
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
        self._wavelength = wavelength
        self._transmittance = transmittance
        self.model = model
        self.filter_set = filter_set
        self.filter_name = filter_name

    @property
    def wavelength(self):
        if self._wavelength is None:
            self._get_filter_from_name()
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        if value is not None:
            self._wavelength = value

    @property
    def transmittance(self):
        if self._transmittance is None:
            self._get_filter_from_name()
        return self._transmittance

    @transmittance.setter
    def transmittance(self, value):
        if value is not None:
            self._transmittance = value

    def _get_filter_from_name(self):
        name = f"{self.filter_set}/{self.filter_name.replace('__', '.')}"
        filt = _filter_generator.reconstruct(name)
        self.wavelength = filt.wavelength
        self.transmittance = filt.transmittance

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

    def __repr__(self):
        if None not in (self.filter_name, self.filter_set):
            return f"<Filter: {self.filter_set}/{self.filter_name}>"
        return "<Filter>"

    def plot(self, ax=None, x_unit=u.um, y_label='Transmittance', **kwargs):

        if ax is None:
            ax = plt.gca()

        label = kwargs.pop('label', self.filter_name)

        with quantity_support():
            ax.plot(
                self.wavelength.to(x_unit), self.transmittance,
                label=label, **kwargs
            )
            ax.legend()

        if y_label is not None:
            ax.set(ylabel=y_label)
