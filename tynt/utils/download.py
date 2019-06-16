from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import download_file
import numpy as np

__all__ = ['DownloadManager']


class DownloadManager(object):
    """
    Manager for downloads from the SVO filter service.

    Examples
    --------
    >>> from tynt import DownloadManager
    >>> # Initialize the filter generator:
    >>> dm = DownloadManager()

    >>> # Included filter sets in download:
    >>> print(dm.include_facilities, dm.include_photsys) # doctest: +SKIP

    >>> # Download all of the links to the filter transmission curves
    >>> dm.download_all_links() # doctest: +SKIP

    >>> # Download all of the filter transmission curve tables
    >>> dm.download_all_tables() # doctest: +SKIP

    >>> # Take the Fourier transform of each transmission curve,
    >>> # output to a FITS BinTable object
    >>> bintable = dm.fft_table() # doctest: +SKIP

    >>> # Write out the BinTable object to a FITS file
    >>> bintable.writeto('fft.fits') # doctest: +SKIP
    """
    include_facilities = ['2MASS', 'SLOAN', 'Kepler', 'TESS', 'HST', 'JWST',
                          'LSST', 'Keck', 'WISE', 'WFIRST', 'Spitzer', 'GAIA']

    include_photsys = ['Bessel', 'Johnson', 'Cousins']

    def __init__(self):
        self.links = None
        self.tables = None

    def download_all_options(self):
        """
        Download all options for instruments, facilities and photometric
        systems from the SVO service.

        This function requires Beautiful Soup as a dependency.

        Returns
        -------
        instruments : list
            List of available instruments
        facilities : list
            List of available facilities
        photometric_system : list
            List of available photometric systems
        """
        from bs4 import BeautifulSoup
        bs = BeautifulSoup(open(download_file('http://svo2.cab.inta-csic.es/'
                                              'theory/fps3/fps.php')).read())
        all_options = [b.attrs['value'] for b in bs.find_all('option')]
        instruments = all_options[1:all_options.index('', 1)]
        facilities = all_options[all_options.index('', 1) + 1:
                                 all_options.index('', 150)]
        photometric_system = all_options[all_options.index('', 150) + 1:]
        return instruments, facilities, photometric_system

    def download_all_links(self, cache=True):
        """
        Download all links to the transmittance curves from the SVO service
        (step 1).

        Parameters
        ----------
        cache : bool (optional)
            Cache the links to your local astropy cache.
        """
        filters = []

        for facility in self.include_facilities:
            url = ('http://svo2.cab.inta-csic.es/theory/fps3/'
                   'fps.php?Facility={0}'
                   .format(facility.replace(" ", "%20")))
            table = Table.read(download_file(url, timeout=300,
                                             cache=cache),
                               format='votable')
            filters.append(table)

        for photsys in self.include_photsys:
            url = ('http://svo2.cab.inta-csic.es/theory/fps3/'
                   'fps.php?PhotSystem={0}'
                   .format(photsys.replace(" ", "%20")))
            table = Table.read(download_file(url, timeout=300,
                                             cache=cache), format='votable')
            filters.append(table)

        self.links = [[i for i in filt['TrasmissionCurve']]
                      for filt in filters]

    def download_all_tables(self, cache=True):
        """
        Download transmittance curve tables (step 2).

        Parameters
        ----------
        cache : bool (optional)
            Cache the links to your local astropy cache.
        """
        if self.links is None:
            raise ValueError("You must run `download_all_links` before "
                             "calling `download_all_tables`.")
        exclude_list = ['Scorpio/Comet.CO+']
        tables = dict()
        for facility_links in self.links:
            for link in facility_links:
                name = link.decode().split('=')[1]
                if name not in tables.keys() and name not in exclude_list:
                    path = download_file(link.decode().replace('+', '&#43;'),
                                         cache=cache)
                    tables[name] = Table.read(path)

        self.tables = tables

    def fft_table(self):
        """
        Generate a FITS BinTable object with the wavelength metadata and
        complex Fourier coefficients representing the filter transmittance
        curves (step 3).

        Returns
        -------
        bt : `~astropy.io.fits.BinTableHDU`
            BinTable object storing complex Fourier coefficients and wavelength
            metadata.
        """
        if self.tables is None:
            raise ValueError("You must run `download_all_tables` before "
                             "calling `fft_table`.")

        rows = dict()
        n_terms = 10

        for k, v in self.tables.items():
            wl, tr = v['Wavelength'], v['Transmission']

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
