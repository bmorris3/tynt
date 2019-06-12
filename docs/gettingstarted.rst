.. include:: references.txt

***************
Getting Started
***************

Approximating a filter transmittance curve
------------------------------------------

First, let's simply plot the approximate transmittance curve SDSS r' filter.
We begin by importing `tynt` and `matplotlib`::

    from tynt import Filter
    import matplotlib.pyplot as plt

We next initialize a `~tynt.Filter` object, which has the useful methods on it::

    f = Filter()

and let's specify the filter that we're observing with::

    identifier = 'SLOAN/SDSS.rprime_filter'

Let's get the approximate transmittance curve as a function of wavelength::

    wavelength, transmittance = f.reconstruct(identifier)

And plot it with `matplotlib`::

    plt.plot(wavelength, transmittance, label=identifier)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Transmittance')
    plt.legend()

.. plot::

    from tynt import Filter
    import matplotlib.pyplot as plt

    f = Filter()

    identifier = 'SLOAN/SDSS.rprime_filter'
    wl, tr = f.reconstruct(identifier)
    plt.plot(wl, tr, label=identifier)
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Transmittance')
    plt.legend()

Plotting all SDSS curves
------------------------

We can grab all of the SDSS prime filters with the following syntax, taking
advantage of the `~tynt.Filter.available_filters()` method, like so::

    import matplotlib.pyplot as plt
    import numpy as np

    from tynt import Filter

    f = Filter()
    filters = [filt for filt in f.available_filters()
               if 'SLOAN/SDSS' in filt and 'prime' in filt]

    fig, ax = plt.subplots(figsize=(10, 5))
    for filt in filters:
        wl, tr = f.reconstruct(filt)
        plt.plot(wl, tr)

        flux_weighted_wl = np.average(wl, weights=tr)

        plt.annotate(filt, xy=(flux_weighted_wl, 0.8 * tr.max()), #filt.split('.')[1]
                     rotation=90, ha='center')
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Transmittance')

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np

    from tynt import Filter

    f = Filter()
    filters = [filt for filt in f.available_filters()
               if 'SLOAN/SDSS' in filt and 'prime' in filt]

    fig, ax = plt.subplots(figsize=(10, 5))
    for filt in filters:
        wl, tr = f.reconstruct(filt)
        plt.plot(wl, tr)

        flux_weighted_wl = np.average(wl, weights=tr)

        plt.annotate(filt, xy=(flux_weighted_wl, 0.8 * tr.max()), #filt.split('.')[1]
                     rotation=90, ha='center')
    plt.xlabel('Wavelength [$\AA$]')
    plt.ylabel('Transmittance')

You can see in the figure above that the Fourier transform approximation does a
rather poor job at the blue-end of the SDSS z' filter.

Comparing the approximation to the true transmittance
-----------------------------------------------------

Finally, let's compare the approximate transmittance curve to the true
transmittance curve, which we'll download from the SVO service::

    import matplotlib.pyplot as plt
    import numpy as np

    from tynt import Filter

    f = Filter()

    filt = 'SLOAN/SDSS.rprime_filter'

    approx_wl, approx_tr = f.reconstruct(filt)
    true_wl, true_tr = f.download_true_transmittance(filt)

    fig, ax = plt.subplots(2, 1, figsize=(4, 8))
    ax[0].plot(true_wl, true_tr, label='True')
    ax[0].plot(approx_wl, approx_tr, label='Approx')
    for axis in ax:
        axis.set_xlabel("Wavelength [$\AA$]")
    ax[0].set_ylabel("Transmittance")
    ax[0].legend()
    ax[1].plot(true_wl, 100*(np.interp(true_wl, approx_wl, approx_tr) - true_tr))
    ax[1].set_ylabel('Error (%)')

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np

    from tynt import Filter

    f = Filter()

    filt = 'SLOAN/SDSS.rprime_filter'

    approx_wl, approx_tr = f.reconstruct(filt)
    true_wl, true_tr = f.download_true_transmittance(filt)

    fig, ax = plt.subplots(2, 1, figsize=(4, 8))
    ax[0].plot(true_wl, true_tr, label='True')
    ax[0].plot(approx_wl, approx_tr, label='Approx')
    for axis in ax:
        axis.set_xlabel("Wavelength [$\AA$]")
    ax[0].set_ylabel("Transmittance")
    ax[0].legend()
    ax[1].plot(true_wl, 100*(np.interp(true_wl, approx_wl, approx_tr) - true_tr))
    ax[1].set_ylabel('Error (%)')

    for axis in ax:
        for s in ['right', 'top']:
            axis.spines[s].set_visible(False)

You can see that the error on the approximate transmittance curve is generally
less than 5% for this filter.