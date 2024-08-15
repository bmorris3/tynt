import numpy as np

from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import download_file

__all__ = ['DownloadManager']

n_terms = 10


class DownloadManager:
    """
    Manager for downloads from the SVO filter service.
    """
    # default facilities included in download:
    include_facilities = ['2MASS', 'SLOAN', 'Kepler', 'TESS', 'HST', 'JWST',
                          'LSST', 'Keck', 'WISE', 'WFIRST', 'Roman', 'Spitzer', 'GAIA',
                          'CHEOPS']

    # default photometric systems included in download:
    include_photsys = ['Bessel', 'Johnson', 'Cousins', 'Stromgren']

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
                                              'theory/fps3/fps.php')).read(), 'xml')
        names = ["INPUT:Instrument", "INPUT:Facility", "INPUT:PhotSystem"]
        instruments, facilities, photometric_system = [
            [
                opt.get("value") for opt in
                bs.findChildren(
                    "PARAM", attrs=dict(name=name)
                )[0].find_all("OPTION")
            ] for name in names
        ]
        return instruments, facilities, photometric_system

    def download_all_links(self, cache=True):
        """
        Download all links to the transmittance curves from the SVO service
        (step 1).

        Parameters
        ----------
        cache : bool
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
        cache : bool
            Cache the links to your local astropy cache.
        """
        if self.links is None:
            raise ValueError("You must run `download_all_links` before "
                             "calling `download_all_tables`.")
        exclude_list = ['Scorpio/Comet.CO+']
        tables = dict()
        for facility_links in self.links:
            for link in facility_links:
                if hasattr(link, 'decode'):
                    link = link.decode()
                name = link.split('=')[1]
                if name not in tables.keys() and name not in exclude_list:
                    if hasattr(link, 'decode'):
                        link = link.decode()
                    path = download_file(link.replace('+', '&#43;'),
                                         cache=cache)
                    tables[name] = Table.read(path)

        self.tables = tables

    def fft_table(self):
        """
        Generate a FITS BinTable object with the wavelength metadata and
        complex Fourier coefficients representing the filter transmittance
        curves (step 3).

        Parameters
        ----------
        n_terms : int
            Number of FFT terms to save

        Returns
        -------
        bt : ~astropy.io.fits.BinTableHDU
            BinTable object storing complex Fourier coefficients and wavelength
            metadata.
        """
        if self.tables is None:
            raise ValueError("You must run `download_all_tables` before "
                             "calling `fft_table`.")

        d = dict()
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
            d[k] = row

        rows = [[id] + d[id] for i, id in enumerate(d.keys())
                if len(d[id]) == n_terms + 4]
        names = (
            "Filter name, n_lambda, lambda_0, delta_lambda, tr_max".split(", ") +
            [f"fft_{i}" for i in range(n_terms)]
        )
        comments = (
            ("SVO FPS filter name, Number of points in transmittance curve, " +
             "Central wavelength (Angstrom), Wavelength spacing (Angstrom), " +
             "Maximum transmittance").split(", ") +
            [f"FFT term: {i}" for i in range(n_terms)]
        )

        filtered_table = Table(rows=rows, names=names)

        bt = fits.BinTableHDU(data=filtered_table)
        for i in range(len(names)):
            bt.header.comments[f"TTYPE{i+1:d}"] = comments[i]

        return bt
