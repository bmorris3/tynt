**********************
Build a custom archive
**********************

Build your own copy of the default tynt filter coefficients
-----------------------------------------------------------

You can build a custom local filter archive by using the :py:class:`~tynt.DownloadManager`
object, like so:

.. code-block:: python

    from tynt import DownloadManager
    # Initialize the filter generator:
    dm = DownloadManager()

    # Included filter sets in download:
    print(dm.include_facilities, dm.include_photsys)

    # Download all of the links to the filter transmission curves
    dm.download_all_links()

    # Download all of the filter transmission curve tables
    dm.download_all_tables()

    # Take the Fourier transform of each transmission curve,
    # output to a FITS BinTable object
    bintable = dm.fft_table()

    # Write out the BinTable object to a FITS file
    bintable.writeto('fft.fits')

Build your own copy of ALL of the SVO FPS filter coefficients
-------------------------------------------------------------

.. note::

    You should be prepared to wait about an hour for the following example to run.
    As of late 2022, there are more than 10,300 filters accessible with this workflow.
    The resulting FITS archive will be ~2.4 MB, and it can zip down to ~1.7 MB.

If you want to include all of the available filters, use the following pattern:

.. code-block:: python

    from tynt import DownloadManager
    # Initialize the filter generator:
    dm = DownloadManager()

    # Find all options available for download (requires Beautiful Soup)
    instruments, facilities, phot_systems = dm.download_all_options()

    dm.include_facilities = facilities
    dm.include_photsys = phot_systems

    # Included filter sets in download:
    print(dm.include_facilities, dm.include_photsys)

    # Download all of the links to the filter transmission curves
    dm.download_all_links()

    # Download all of the filter transmission curve tables
    dm.download_all_tables()

    # Take the Fourier transform of each transmission curve,
    # output to a FITS BinTable object
    bintable = dm.fft_table()

    # Write out the BinTable object to a FITS file
    bintable.writeto('fft.fits')

Now when you initialize the :py:class:`~tynt.FilterGenerator` object, you can supply
it with the path to your newly created ``fft.fits`` file.