import numpy as np
from astropy.tests.helper import remote_data, assert_quantity_allclose
import astropy.units as u

from ..core import FilterGenerator

f = FilterGenerator()


@remote_data
def test_download_sdss_r():

    identifier = 'SLOAN/SDSS.rprime_filter'

    filt = f.reconstruct(identifier)
    approx_wl, approx_tr = filt.wavelength, filt.transmittance
    true_filt = f.download_true_transmittance(identifier)
    true_wl, true_tr = true_filt.wavelength, true_filt.transmittance

    np.testing.assert_allclose(np.interp(true_wl, approx_wl, approx_tr),
                               true_tr, atol=0.5)


@remote_data
def test_download_sdss_g():

    identifier = 'SLOAN/SDSS.gprime_filter'

    filt = f.reconstruct(identifier)
    approx_wl, approx_tr = filt.wavelength, filt.transmittance
    true_filt = f.download_true_transmittance(identifier)
    true_wl, true_tr = true_filt.wavelength, true_filt.transmittance

    np.testing.assert_allclose(np.interp(true_wl, approx_wl, approx_tr),
                               true_tr, atol=0.07)


@remote_data
def test_download_sdss_i():

    identifier = 'SLOAN/SDSS.iprime_filter'

    filt = f.reconstruct(identifier)
    approx_wl, approx_tr = filt.wavelength, filt.transmittance
    true_filt = f.download_true_transmittance(identifier)
    true_wl, true_tr = true_filt.wavelength, true_filt.transmittance

    np.testing.assert_allclose(np.interp(true_wl, approx_wl, approx_tr),
                               true_tr, atol=0.05)


def test_sdss_lambda_eff_w_eff():

    sdss_filters = ['SLOAN/SDSS.u', 'SLOAN/SDSS.g',
                    'SLOAN/SDSS.r', 'SLOAN/SDSS.i',
                    'SLOAN/SDSS.z']

    # Answers from:
    # http://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=SLOAN
    lambda_mean_true = u.Quantity([3561.8, 4718.9, 6185.2, 7499.7, 8961.5],
                                  u.Angstrom)

    w_eff_true = u.Quantity([558.4, 1158.4, 1111.2, 1044.6, 1124.6],
                             u.Angstrom)
    w_eff_approx = []
    lambda_mean_approx = []

    for identifier in sdss_filters:
        filt = f.reconstruct(identifier)

        lambda_bar_approx = (np.trapz(filt.transmittance * filt.wavelength,
                                      filt.wavelength) /
                             np.trapz(filt.transmittance, filt.wavelength))

        width_approx = (np.trapz(filt.transmittance, filt.wavelength) /
                        filt.transmittance.max())

        w_eff_approx.append(width_approx)
        lambda_mean_approx.append(lambda_bar_approx)

    assert_quantity_allclose(w_eff_approx, w_eff_true,
                             rtol=0.05)
    assert_quantity_allclose(lambda_mean_approx, lambda_mean_true, rtol=0.01)
