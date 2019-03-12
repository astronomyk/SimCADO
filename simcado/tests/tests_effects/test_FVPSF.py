# 1. FVPSF should return the PSF for a position in the FOV and a given lambda
# 2. should throw errors when:
#   - file doesn't exist
#   - file doesn't have a image or table in the 1st extension
# 3. should have attributes:
#   - lam_bin_centers : pulled from the header of the hduNs
#   - layer_map : an ImageHDU with a map of where each layer is valid
#   - layer_table : The BinTableHDU if available
#   - _make_layer_map : makes a layer_map from a BinTableHDU is needed
#   - mask(wave, pos) : returns an ImageHDU with WCS with a mask of the valid
#                       region for a given wavelength and position in the FOV
#   - shape : (N_EXT, N_LAYER, NAXIS2, NAXIS1)
#   - nearest(wave, pos=None, hdu=False) : should return the array for the
#                                          given position and wavelength
#   - defaults : a dictionary with {"wave": , pos: (,)} so that .array can be
#                used for backwards compatibility
#   - set_defaults(wave, pos)
#   - array : returns an array for the defaults values of wave and pos
#   - psf : returns self, as array returns a layer based on defaults

import os
from copy import deepcopy
import pytest
from pytest import approx

import numpy as np
from astropy import units as u
from astropy.io import fits

import simcado as sim
from simcado.optics.fov import FieldOfView
from simcado.optics import image_plane_utils as imp_utils
from simcado.optics.optical_train import OpticalTrain
from simcado.utils import find_file
from simcado.commands.user_commands2 import UserCommands
from simcado.optics.effects.psfs import FieldVaryingPSF
from simcado.optics.effects import psfs

from simcado.tests.mocks.py_objects.source_objects import _image_source, \
    _single_table_source
from simcado.tests.mocks.py_objects.psf_objects import _basic_circular_fvpsf

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
PLOTS = False

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))
sim.rc.__search_path__ += [FILES_PATH, YAMLS_PATH]


def _centre_fov():
    n = 55
    xsky = np.array([-n, n]) * u.arcsec.to(u.deg)
    ysky = np.array([-n, n]) * u.arcsec.to(u.deg)
    sky_hdr = imp_utils.header_from_list_of_xy(xsky, ysky, 1/3600.)
    imp_hdr = imp_utils.header_from_list_of_xy([-n, n], [-n, n], 1, "D")
    imp_hdr.update(sky_hdr)
    return FieldOfView(imp_hdr, waverange=[1.0, 2.0]*u.um)


@pytest.fixture(scope="function")
def centre_fov():
    return _centre_fov()


@pytest.fixture(scope="function")
def basic_circular_fvpsf():
    return _basic_circular_fvpsf


class TestInit:
    def test_errors_when_initialised_with_nothing(self):
        with pytest.raises(ValueError):
            isinstance(FieldVaryingPSF(), FieldVaryingPSF)

    def test_initialised_when_passed_fits_filename(self):
        fvpsf = FieldVaryingPSF(filename="test_FVPSF.fits")
        assert isinstance(fvpsf, FieldVaryingPSF)


class TestStrehlImageHDU:
    def test_returns_the_correct_hdu(self):
        fvpsf = FieldVaryingPSF(filename="test_FVPSF.fits")
        strehl_imhdu = fvpsf.strehl_imagehdu
        assert isinstance(strehl_imhdu, fits.ImageHDU)


@pytest.mark.usefixtures("centre_fov")
class TestGetKernel:
    def test_returns_array_with_single_kernel_from_fov(self, centre_fov):
        fvpsf = FieldVaryingPSF(filename="test_FVPSF.fits")
        kernels = fvpsf.get_kernel(centre_fov)
        assert np.all(kernels[0][0] == fvpsf._file[2].data[4])
        assert kernels[0][1] is None

    def test_returns_four_arrays_when_fov_on_intersection(self, centre_fov):
        centre_fov.hdu.header["CRVAL1"] -= 15/3600.
        centre_fov.hdu.header["CRVAL2"] -= 15/3600.
        fvpsf = FieldVaryingPSF(filename="test_FVPSF.fits")
        kernels = fvpsf.get_kernel(centre_fov)
        assert len(kernels) == 4

        if PLOTS:
            for ii, kernel_mask in enumerate(kernels):
                plt.subplot(2, 4, ii+1)
                plt.imshow(kernel_mask[0].T, origin="lower")
                plt.subplot(2, 4, ii+5)
                plt.imshow(kernel_mask[1].T, origin="lower")

            plt.show()


@pytest.mark.usefixtures("centre_fov", "basic_circular_fvpsf")
class TestApplyTo:
    def test_convolution_with_delta_for_central_region(self, centre_fov):
        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]
        centre_fov.hdu.data = np.zeros((nax1, nax2))
        centre_fov.hdu.data[::3, ::3] = 1

        fvpsf = FieldVaryingPSF(filename="test_FVPSF.fits")
        fov_back = fvpsf.apply_to(centre_fov)

        if PLOTS:
            plt.imshow(fov_back.hdu.data.T, origin="lower")
            plt.show()

    def test_convolution_with_fvpsfs_for_shifted_region(self, centre_fov):
        # centre_fov.hdu.header["CRVAL1"] -= 15/3600.
        # centre_fov.hdu.header["CRVAL2"] -= 15/3600.

        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]
        centre_fov.hdu.data = np.zeros((nax1, nax2))
        centre_fov.hdu.data[::5, ::5] = 1
        centre_fov.fields = [1]

        fvpsf = FieldVaryingPSF(filename="test_FVPSF.fits")
        fov_back = fvpsf.apply_to(centre_fov)

        if PLOTS:
            plt.imshow(fov_back.hdu.data.T, origin="lower")
            plt.show()

    def test_circular_fvpsf(self, centre_fov, basic_circular_fvpsf):
        # centre_fov.hdu.header["CRVAL1"] -= 15/3600.
        # centre_fov.hdu.header["CRVAL2"] -= 15/3600.

        nax1, nax2 = centre_fov.header["NAXIS1"], centre_fov.header["NAXIS2"]
        centre_fov.hdu.data = np.zeros((nax1, nax2))
        x, y = np.random.randint(0, 100, (2, 150))
        centre_fov.hdu.data[x, y] = 1
        # centre_fov.hdu.data[::11, ::11] = 1
        centre_fov.fields = [1]

        fvpsf = FieldVaryingPSF(filename="test_circular_fvpsf.fits")
        fov_back = fvpsf.apply_to(centre_fov)

        if not PLOTS:
            plt.imshow(fov_back.hdu.data.T, origin="lower", vmax=0.1)
            plt.show()


@pytest.mark.usefixtures("centre_fov")
class TestFunctionGetStrehlCutout:
    @pytest.mark.parametrize("scale", [0.2, 0.5, 1, 2])
    def test_returns_correct_section_of_strehl_map(self, centre_fov, scale):
        centre_fov.hdu.header["CDELT1"] *= scale
        centre_fov.hdu.header["CDELT2"] *= scale

        fvpsf = FieldVaryingPSF(filename="test_FVPSF.fits")
        strehl_hdu = psfs.get_strehl_cutout(centre_fov.header,
                                            fvpsf.strehl_imagehdu)

        #assert all(np.unique(strehl_hdu.data).astype(int) == [0, 1, 3, 4])

        if not PLOTS:
            plt.imshow(strehl_hdu.data.T, origin="lower")
            plt.colorbar()
            plt.show()
