from copy import deepcopy
import numpy as np
import pytest
from pytest import approx

from astropy import units as u

from simcado.optics.fov import FieldOfView
from simcado.optics.image_plane import ImagePlane

from simcado.tests.mocks.py_objects import source_objects as src
from simcado.tests.mocks.py_objects import header_objects as hdrs

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False


@pytest.fixture(scope="function")
def comb_src():
    return src._combined_source()


@pytest.fixture(scope="function")
def fov_hdr():
    return hdrs._fov_header()


@pytest.fixture(scope="function")
def implane_hdr():
    return hdrs._implane_header()


@pytest.mark.usefixtures("comb_src", "fov_hdr", "implane_hdr")
class TestCanExtractSourceAndPlaceOnImagePlane:
    def test_can_extract_the_source_in_a_fov(self, fov_hdr, comb_src,
                                             implane_hdr):
        fov_hdr["CRVAL1D"] += 20
        fov_hdr["CRVAL2D"] += 20

        fov = FieldOfView(fov_hdr, waverange=[0.5, 2.5]*u.um)
        imp = ImagePlane(implane_hdr)

        fov.extract_from(comb_src)
        imp.add(fov.fields, wcs_suffix="D")

        ipt = np.sum(fov.fields[0]["flux"]) + np.sum(fov.fields[1].data)
        opt = np.sum(imp.image)

        assert ipt == approx(opt)

        if PLOTS:
            plt.subplot(121)
            plt.imshow(fov.image.T, origin="lower", norm=LogNorm())

            plt.subplot(122)
            plt.imshow(imp.image.T, origin="lower", norm=LogNorm())
            plt.show()




