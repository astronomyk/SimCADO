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
# Helper functions

def mock_dir():
    cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
    rel_dirname = "mocks"
    return os.path.abspath(os.path.join(cur_dirname, rel_dirname))

MOCK_DIR = mock_dir()


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
    This functions measures the counts within circular apertures

    Parameters
    ----------
    image: a image

    Returns
    -------
    counts
    """
    (X, Y) = (512, 512)
    reff = 20
    indices = np.indices(np.shape(image))
    distances = np.sqrt((X - indices[0]) ** 2 + (Y - indices[1]) ** 2)
    aperture_mask = np.ma.masked_less_equal(distances, 6 * reff)
    sky_mask = np.ma.masked_greater_equal(distances, 12 * reff)
    aperture_area = np.size(image) - aperture_mask.count()
    sky_area = np.size(image) - sky_mask.count()
    aperture_flux = aperture_mask.mask * image
    sky_flux = sky_mask.mask * image
    counts = np.sum(aperture_flux) - np.sum(sky_flux) * (aperture_area / sky_area)
    return counts

# End helper functions


def test_create_image_scaled_by_factor(factor=10):
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


def test_photometry(sky_level=100,factor=10):
    """
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


def test_source_from_image(factor1=10, factor2=20):
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
    galaxy_src1 = sim.source.source_from_image(image1, lam, spec, plate_scale=0.004, pix_res=0.004,
                                              flux_threshold=0, conserve_flux=True)
    galaxy_src2 = sim.source.source_from_image(image2, lam, spec, plate_scale=0.004, pix_res=0.004,
                                              flux_threshold=0, conserve_flux=True)
    sim_img1 = sim.run(galaxy_src1, OBS_NDIT=1, OBS_EXPTIME=100, SIM_DETECTOR_PIX_SCALE=0.004)
    sim_img2 = sim.run(galaxy_src2, OBS_NDIT=1, OBS_EXPTIME=100, SIM_DETECTOR_PIX_SCALE=0.004)
    counts1 = photometry(sim_img1[0].data)
    counts2 = photometry(sim_img2[0].data)
    assert np.abs(counts1/counts2 - 1) < 0.1


def test_source_elliptical():
    """
    Test whether source.elliptical produces consistent results in comparison with source from image

    Returns
    -------

    """
    filter_file = os.path.join(MOCK_DIR, 'TC_filter_K.dat')
    image1 = create_image_scaled_by_factor(1)
    lam, spec = sim.source.SED("spiral", filter_name=filter_file, magnitude=15)
    galaxy_src1 = sim.source.source_from_image(image1, lam, spec, plate_scale=0.004, pix_res=0.004,
                                              flux_threshold=0, conserve_flux=True)
    galaxy_src2 = sim.source.elliptical(20 * 0.004, 0.004, magnitude=15, n=1,
                                        filter_name=filter_file, normalization="total",
                                        spectrum="spiral", ellipticity=0, angle=0)
    sim_img1 = sim.run(galaxy_src1, OBS_NDIT=1, OBS_EXPTIME=100, SIM_DETECTOR_PIX_SCALE=0.004)
    sim_img2 = sim.run(galaxy_src2, OBS_NDIT=1, OBS_EXPTIME=100, SIM_DETECTOR_PIX_SCALE=0.004)
    counts1 = photometry(sim_img1[0].data)
    counts2 = photometry(sim_img2[0].data)
    assert np.abs(counts1 / counts2 - 1) < 0.1


#
#
