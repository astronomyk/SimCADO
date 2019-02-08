# Integration test functionality for simcado.source.source_from_image()
#
# 1) Two images are created scaled by arbitrary factors
# 2) source_from_image should scale both of them at the same level
#   according to the information provided by source.SED
# 3) Perform aperture photometry to check that number counts are consistent
#
# Source is a n=1 (exponential profile) with Reff = 20 pixels and zero ellipticity

import numpy as np
import simcado as sim
import os
import inspect
import pytest
from simcado import __data_dir__

# Helper functions ---

def mock_dir():
    cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
    rel_dirname = "mocks"
    return os.path.abspath(os.path.join(cur_dirname, rel_dirname))

MOCK_DIR = mock_dir()

sim.get_extras()
cmds = sim.UserCommands(os.path.join(MOCK_DIR, "default.conf"))
cmds["SIM_DATA_DIR"] = __data_dir__
print(cmds["SIM_DATA_DIR"])
print(cmds["ATMO_EC"])
# Update some keywords by hand
cmds["INST_ENTR_WINDOW_TC"] = os.path.join(MOCK_DIR, "TC_window.dat")
cmds["INST_DICHROIC_TC"] = os.path.join(MOCK_DIR, "TC_dichroic.dat")
cmds["INST_PUPIL_TC"] = os.path.join(MOCK_DIR, "TC_pupil.dat")
cmds["INST_ADC_TC"] = os.path.join(MOCK_DIR, "TC_ADC.dat")


sim.get_extras()


def create_image_scaled_by_factor(factor=1):
    """
    This creates a simple source with a SERSIC profile  multiplied by a factor
    to simulate a random input source


    NOTE: n=1 to make the light profile decline fast and so the rough photometry function works

    Parameters
    ----------
    factor = a multiplicative number (sersic profile images are normalized to 1)

    Returns
    -------
    image * factor
    """
    image = sim.source.sersic_profile(r_eff=20, n=1, ellipticity=0., angle=0,
                                      normalization="total", width=1024, height=1024, x_offset=0, y_offset=0,
                                      oversample=1)
    image = image * factor
    return image


def photometry(image):
    """
    Very simple function to return the counts of a single object
    it sums the counts and substract a sky level

    Parameters
    ----------
    image: a image

    Returns
    -------
    counts
    """

    npix = np.size(image)
    sky = np.median(image[0:100,:])
    total = np.sum(image)
    obj = total - sky * npix
    return obj

# End of helper functions ---


@pytest.mark.parametrize("factor", np.arange(20))
def test_create_image_scaled_by_factor(factor):
    """
    Test if the image is created as expected
    Parameters
    ----------
    factor: a multiplicative factor to the image

    Returns
    -------

    """
    image = create_image_scaled_by_factor(factor)
    counts = np.sum(image)
    assert np.abs(counts - factor) < 0.01


@pytest.mark.parametrize("sky_level", np.linspace(1, 25, 5)**3)
@pytest.mark.parametrize("factor", np.linspace(1, 25, 5))
def test_photometry(sky_level, factor):
    """re
    Test if photometry function is performing as expected

    Parameters
    ----------
    sky_level: an user chosen background level
    factor: multiplicative factor passed to create_image_scaled_by_factor

    Returns
    -------

    """
    image = create_image_scaled_by_factor(factor) + sky_level
    counts = photometry(image)
    assert np.abs(counts - factor) < 0.01


@pytest.mark.parametrize("factor1", np.linspace(1, 25, 5))
@pytest.mark.parametrize("factor2", np.linspace(1, 25, 5))
def test_source_from_image(factor1, factor2):
    """
    test source from image

    Two images are created scaled by different factors. The images are supposed to be scaled
    to the specified magnitude by source_from_image() according to the information provided by source.SED()

    Test the counts are consistent after a simulation.

    Parameters
    ----------
    factor1: factor to create image1
    factor2: factor to create image2

    Returns
    -------

    """

    filter_file = os.path.join(MOCK_DIR, 'TC_filter_K.dat')
    image1 = create_image_scaled_by_factor(factor1)
    image2 = create_image_scaled_by_factor(factor2)
    lam, spec = sim.source.SED("spiral", filter_name=filter_file, magnitude=15)
    galaxy_src1 = sim.source.source_from_image(image1, lam, spec, plate_scale=0.004,
                                               pix_res=0.004, flux_threshold=0, conserve_flux=True)
    galaxy_src2 = sim.source.source_from_image(image2, lam, spec, plate_scale=0.004,
                                               pix_res=0.004, flux_threshold=0, conserve_flux=True)
    sim_img1 = sim.run(galaxy_src1, OBS_NDIT=1, OBS_EXPTIME=300, SIM_DETECTOR_PIX_SCALE=0.004, cmds=cmds)
    sim_img2 = sim.run(galaxy_src2, OBS_NDIT=1, OBS_EXPTIME=300, SIM_DETECTOR_PIX_SCALE=0.004, cmds=cmds)
    counts1 = photometry(sim_img1[0].data)
    counts2 = photometry(sim_img2[0].data)
    assert np.abs(counts1 / counts2 - 1) < 0.1


@pytest.mark.parametrize("mag", np.linspace(9, 18, 10))
def test_source_elliptical(mag):
    """
    Test whether source.elliptical produces consistent results in comparison with source from image

    Returns
    -------

    """
    filter_file = os.path.join(MOCK_DIR, 'TC_filter_K.dat')
    image1 = create_image_scaled_by_factor(1)
    lam, spec = sim.source.SED("spiral", filter_name=filter_file, magnitude=mag)
    galaxy_src1 = sim.source.source_from_image(image1, lam, spec, plate_scale=0.004,
                                               pix_res=0.004, flux_threshold=0, conserve_flux=True)
    galaxy_src2 = sim.source.elliptical(20 * 0.004, 0.004, magnitude=mag, n=1,
                                        filter_name=filter_file, normalization="total",
                                        spectrum="spiral", ellipticity=0, angle=0)
    sim_img1 = sim.run(galaxy_src1, OBS_NDIT=1, OBS_EXPTIME=300, SIM_DETECTOR_PIX_SCALE=0.004, cmds=cmds)
    sim_img2 = sim.run(galaxy_src2, OBS_NDIT=1, OBS_EXPTIME=300, SIM_DETECTOR_PIX_SCALE=0.004, cmds=cmds)
    counts1 = photometry(sim_img1[0].data)
    counts2 = photometry(sim_img2[0].data)
    print(mag)
    assert np.abs(counts1 / counts2 - 1) < 0.1



#
#
