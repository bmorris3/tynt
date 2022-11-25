====
tynt
====

``tynt`` is a super lightweight package containing approximate transmittance
curves for more than one thousand astronomical filters, weighing in at just
186 KB. Depends only on numpy and astropy.

Contents
--------

.. toctree::
  :maxdepth: 2

  tynt/installation
  tynt/gettingstarted
  tynt/customarchive
  tynt/filters.rst
  tynt/api.rst

Simple example
^^^^^^^^^^^^^^

Let's plot the transmittance curve of the SDSS r' filter:

.. plot::
    :include-source:

    from tynt import FilterGenerator

    f = FilterGenerator()
    filt = f.reconstruct('SLOAN/SDSS.rprime_filter')

    import matplotlib.pyplot as plt
    plt.plot(filt.wavelength, filt.transmittance)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Approx. Transmittance')

Links
^^^^^

* `Source code (GitHub) <https://github.com/bmorris3/tynt>`_
* `Docs <https://tynt.readthedocs.io/>`_
* `Issues <https://github.com/bmorris3/tynt/issues>`_

References
^^^^^^^^^^

This research has made use of the SVO Filter Profile Service
(http://svo2.cab.inta-csic.es/theory/fps/) supported from the Spanish MINECO
through grant AYA2017-84089

* `The SVO Filter Profile Service. Rodrigo, C., Solano, E., Bayo, A. <http://ivoa.net/documents/Notes/SVOFPS/index.html>`_
* `The Filter Profile Service Access Protocol. Rodrigo, C., Solano, E. <http://ivoa.net/documents/Notes/SVOFPSDAL/index.html>`_
