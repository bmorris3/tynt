***************
Getting Started
***************

Approximating a filter transmittance curve
------------------------------------------

First, let's simply plot the approximate transmittance curve SDSS r' filter.
We begin by importing ``tynt`` and ``matplotlib``:

.. code-block:: python

    from tynt import FilterGenerator
    import matplotlib.pyplot as plt

We next initialize a :py:class:`~tynt.FilterGenerator` object, which has the useful methods on it:

.. code-block:: python

    f = FilterGenerator()

and let's specify the filter that we're observing with:

.. code-block:: python

    identifier = 'SLOAN/SDSS.rprime_filter'

Let's get the approximate transmittance curve as a function of wavelength:

.. code-block:: python

    filt = f.reconstruct(identifier)

And plot it with ``matplotlib``:

.. code-block:: python

    import matplotlib.pyplot as plt
    plt.plot(filt.wavelength.value, filt.transmittance, label=identifier)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Transmittance')
    plt.legend()

.. plot::

    from tynt import FilterGenerator
    import matplotlib.pyplot as plt

    f = FilterGenerator()

    identifier = 'SLOAN/SDSS.rprime_filter'
    filt = f.reconstruct(identifier)

    plt.plot(filt.wavelength.value, filt.transmittance, label=identifier)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Transmittance')
    plt.legend()

Plotting all SDSS curves
------------------------

We can grab all of the SDSS prime filters with the following syntax, taking
advantage of the :py:meth:`~tynt.FilterGenerator.available_filters` method, like so:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import numpy as np

    from tynt import FilterGenerator

    f = FilterGenerator()
    filters = [filt for filt in f.available_filters()
               if 'SLOAN/SDSS' in filt and 'prime' in filt]

    fig, ax = plt.subplots(figsize=(10, 5))

    for filt in filters:
        sdss_filter = f.reconstruct(filt)
        plt.plot(sdss_filter.wavelength.value, sdss_filter.transmittance)

        flux_weighted_wl = np.average(sdss_filter.wavelength.value,
                                      weights=sdss_filter.transmittance)

        plt.annotate(filt, xy=(flux_weighted_wl,
                               0.8 * sdss_filter.transmittance.max()),
                     rotation=90, va='top')

    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Transmittance')


You can see in the figure above that the Fourier transform approximation does a
rather poor job at the blue-end of the SDSS z' filter.

Comparing the approximation to the true transmittance
-----------------------------------------------------

Finally, let's compare the approximate transmittance curve to the true
transmittance curve, which we'll download from the SVO service:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    import numpy as np

    from tynt import FilterGenerator

    f = FilterGenerator()

    filt = 'SLOAN/SDSS.rprime_filter'

    filt_approx = f.reconstruct(filt)
    filt_true = f.download_true_transmittance(filt)

    fig, ax = plt.subplots(2, 1, figsize=(4, 8))
    ax[0].plot(filt_true.wavelength.value, filt_true.transmittance, label='True')
    ax[0].plot(filt_approx.wavelength.value, filt_approx.transmittance, label='Approx')
    for axis in ax:
        axis.set_xlabel("Wavelength [$\AA$]")
    ax[0].set_ylabel("Transmittance")
    ax[0].legend()

    difference = 100*(np.interp(filt_true.wavelength,
                                filt_approx.wavelength,
                                filt_approx.transmittance) - filt_true.transmittance)

    ax[1].plot(filt_true.wavelength.value, difference)
    ax[1].set_ylabel('Error (%)')

    for axis in ax:
        for s in ['right', 'top']:
            axis.spines[s].set_visible(False)

You can see that the error on the approximate transmittance curve is generally
less than 5% for this filter.

Constructing an astropy model transmittance curve
-------------------------------------------------

In some instances it may be useful to represent the transmittance curve
analytically with an astropy model. You can get a custom astropy model
like so:

.. plot::
    :include-source:

    from tynt import FilterGenerator
    import matplotlib.pyplot as plt

    f = FilterGenerator()

    identifier = 'SLOAN/SDSS.rprime_filter'
    filt = f.reconstruct(identifier, model=True)

    plt.plot(filt.wavelength.value, filt.model(filt.wavelength.value))
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Transmittance')

Getting transmittance curves not included by default
----------------------------------------------------

A list of the filters available in ``tynt`` is documented in :doc:`filters`.
You also have access to *all* of the filters stored in the
`SVO Filter Profile Service <http://svo2.cab.inta-csic.es/theory/fps/>`_ if you
have internet access via the :py:meth:`~tynt.FilterGenerator.download_true_transmittance`
method:

.. plot::
    :include-source:

    from tynt import FilterGenerator
    import matplotlib.pyplot as plt

    f = FilterGenerator()

    identifier_b = 'TYCHO/TYCHO.B'
    identifier_v = 'TYCHO/TYCHO.V'
    filt_b = f.download_true_transmittance(identifier_b)
    filt_v = f.download_true_transmittance(identifier_v)

    plt.plot(filt_b.wavelength.value, filt_b.transmittance)
    plt.plot(filt_v.wavelength.value, filt_v.transmittance)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Transmittance')

To make a local archive of *all* filters available via SVO, see :doc:`customarchive`.
