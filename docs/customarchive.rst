.. include:: references.txt

**********************
Build a custom archive
**********************

You can build a custom local filter archive by using the `~tynt.DownloadManager`
object, like so::

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

If you want to include all of the available filters, use the following pattern::

    >>> from tynt import DownloadManager
    >>> # Initialize the filter generator:
    >>> dm = DownloadManager()

    >>> # Find all options available for download (requires Beautiful Soup)
    >>> instruments, facilities, phot_systems = dm.download_all_options() # doctest: +SKIP

    >>> dm.include_facilities = facilities # doctest: +SKIP
    >>> dm.include_photsys = phot_systems # doctest: +SKIP

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

Now when you initialize the `~tynt.FilterGenerator` object, you can supply
it with the path to your newly created `fft.fits` file.