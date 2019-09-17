"""
Tests for:

simulation.zeropoint
"""

import numpy as np
import simcado
import pytest


@pytest.mark.parametrize("filter_name", ["TC_filter_Ks.dat"])
@pytest.mark.parametrize("mag", np.linspace(10, 20, 3))
@pytest.mark.parametrize("pixel_size", [0.004, 0.0015])
@pytest.mark.parametrize("fwhm", np.linspace(3, 13, 3))
def test_zeropoint(filter_name, mag, pixel_size, fwhm):
    """
    Modified version of source.zeropoint:
    Test that default parameters do not make a difference
    in the returned zeropoint value

    Returns the zero point magnitude for a SimCADO filter

    This is an end-to-end simulation which aims to take into account all transmission effects
    incorporated in a SimCADO simulation.

    The returned zeropoint is for an exposure of 1s. Therefore, magnitudes from measured fluxes in simulated
    images should be calculated as following

    mag = -2.5*np.log10(counts/texp) + zp

    where counts are the background subtracted counts, texp is the exposure time and zp is the zeropoint for
    the filter in question, calculated here.

    Parameters
    ----------
    filter_name: A SimCADO filter

    Returns
    -------
    zp: the zeropoint magnitude


    """
    msg = "Calculating zeropoint for filter:", filter_name
    print(msg)
    input_mag = mag
    pixel_size = pixel_size
    fwhm = fwhm
    seeing = fwhm * pixel_size
    texp = 1
    star = simcado.source.star(spec_type="A0V", mag=input_mag, filter_name=filter_name, x=0, y=0)
    hdu = simcado.run(star, OBS_DIT=texp,
                      INST_FILTER_TC=filter_name, SIM_DETECTOR_PIX_SCALE=pixel_size,
                      SCOPE_PSF_FILE=None, OBS_SEEING=seeing,
                      FPA_LINEARITY_CURVE=None, FPA_CHIP_LAYOUT="small", FPA_USE_NOISE="no",
                      ATMO_USE_ATMO_BG="no", SCOPE_USE_MIRROR_BG="no", INST_USE_AO_MIRROR_BG="no")
    image = hdu[0].data
    sky_level = np.median(image[0:200, :])
    x1, x2 = 512 - 12 * int(fwhm), 512 + 12 * int(fwhm)
    y1, y2 = 512 - 12 * int(fwhm), 512 + 12 * int(fwhm)
    cut_image = image[x1:x2, y1:y2] - sky_level
    counts = np.sum(cut_image)
    zp = 2.5 * np.log10(counts) + input_mag
    zp = np.round(zp, 3)

    assert np.isclose(zp, 29.491, 1e-2)  # Comparing to tabulated ZP


@pytest.mark.parametrize("filter_name", ["TC_filter_J.dat", "TC_filter_Ks.dat"])
@pytest.mark.parametrize("mag_in", [12, 17, 22])
@pytest.mark.parametrize("texp", [10, 100, 1000])
def test_input_magnitudes(filter_name, mag_in, texp):
    """
    Test that input and output magnitudes are the same

    """
    star = simcado.source.star(spec_type="A0V", mag=mag_in,
                               filter_name=filter_name, x=0, y=0)
    pixel_size = 0.004
    fwhm = 3.6
    seeing = fwhm * pixel_size
    hdu = simcado.run(star, OBS_DIT=texp,
                      INST_FILTER_TC=filter_name, SIM_DETECTOR_PIX_SCALE=pixel_size,
                      SCOPE_PSF_FILE=None, OBS_SEEING=seeing,
                      FPA_LINEARITY_CURVE=None, FPA_CHIP_LAYOUT="small", FPA_USE_NOISE="no",
                      ATMO_USE_ATMO_BG="no", SCOPE_USE_MIRROR_BG="no", INST_USE_AO_MIRROR_BG="no")
    image = hdu[0].data
    sky_level = np.median(image[0:200, :])
    x1, x2 = 512 - 12 * int(fwhm), 512 + 12 * int(fwhm)
    y1, y2 = 512 - 12 * int(fwhm), 512 + 12 * int(fwhm)
    cut_image = image[x1:x2, y1:y2] - sky_level
    counts = np.sum(cut_image)

    mag_out = -2.5*np.log10(counts/texp) + simcado.simulation.zeropoint(filter_name)

    assert np.isclose(mag_in, mag_out, 1e-2)



