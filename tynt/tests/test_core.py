from astropy.tests.helper import remote_data
import numpy as np

from ..core import FilterGenerator

f = FilterGenerator()


@remote_data
def test_sdss_r():

    identifier = 'SLOAN/SDSS.rprime_filter'

    filt = f.reconstruct(identifier)
    approx_wl, approx_tr = filt.wavelength, filt.transmittance
    true_filt = f.download_true_transmittance(identifier)
    true_wl, true_tr = true_filt.wavelength, true_filt.transmittance

    np.testing.assert_allclose(np.interp(true_wl, approx_wl, approx_tr),
                               true_tr, atol=0.5)


@remote_data
def test_sdss_g():

    identifier = 'SLOAN/SDSS.gprime_filter'

    filt = f.reconstruct(identifier)
    approx_wl, approx_tr = filt.wavelength, filt.transmittance
    true_filt = f.download_true_transmittance(identifier)
    true_wl, true_tr = true_filt.wavelength, true_filt.transmittance

    np.testing.assert_allclose(np.interp(true_wl, approx_wl, approx_tr),
                               true_tr, atol=0.07)


@remote_data
def test_sdss_i():

    identifier = 'SLOAN/SDSS.iprime_filter'

    filt = f.reconstruct(identifier)
    approx_wl, approx_tr = filt.wavelength, filt.transmittance
    true_filt = f.download_true_transmittance(identifier)
    true_wl, true_tr = true_filt.wavelength, true_filt.transmittance

    np.testing.assert_allclose(np.interp(true_wl, approx_wl, approx_tr),
                               true_tr, atol=0.05)
