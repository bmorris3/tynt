from astropy.tests.helper import remote_data
import numpy as np

from ..core import Filter

f = Filter()


@remote_data
def test_sdss_r():

    filt = 'SLOAN/SDSS.rprime_filter'

    approx_wl, approx_tr = f.reconstruct(filt)
    true_wl, true_tr = f.download_true_transmittance(filt)

    np.testing.assert_allclose(np.interp(true_wl, approx_wl, approx_tr),
                               true_tr, atol=0.5)


@remote_data
def test_sdss_g():

    filt = 'SLOAN/SDSS.gprime_filter'

    approx_wl, approx_tr = f.reconstruct(filt)
    true_wl, true_tr = f.download_true_transmittance(filt)

    np.testing.assert_allclose(np.interp(true_wl, approx_wl, approx_tr),
                               true_tr, atol=0.07)


@remote_data
def test_sdss_i():

    filt = 'SLOAN/SDSS.iprime_filter'

    approx_wl, approx_tr = f.reconstruct(filt)
    true_wl, true_tr = f.download_true_transmittance(filt)

    np.testing.assert_allclose(np.interp(true_wl, approx_wl, approx_tr),
                               true_tr, atol=0.05)
