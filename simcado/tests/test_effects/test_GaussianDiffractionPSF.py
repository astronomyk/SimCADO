import pytest

import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy import units as u

from simcado.optics.effects import gaussian_diffraction_psf as gdf
from simcado.optics.fov import FieldOfView
from simcado.source.source2 import Source

from simcado.tests.mocks.objects.source_objects import _image_source
from simcado.tests.mocks.objects.header_objects import _basic_fov_header


def _make_basic_fov():
    src = _image_source()
    fov = FieldOfView(_basic_fov_header(), waverange=[1, 2]*u.um)
    fov.extract_from(src)

    return fov


class TestInit:
    def test_does_not_initialise_with_nothing(self):
        with pytest.raises(TypeError):
            gdf.GaussianDiffractionPSF()

    def test_initialises_with_only_diameter(self):
        eff = gdf.GaussianDiffractionPSF(1)
        assert isinstance(eff, gdf.GaussianDiffractionPSF)

    def test_initialised_with_other_keywords(self):
        eff = gdf.GaussianDiffractionPSF(1, sub_pixel=False)
        assert eff.meta["sub_pixel"] == False


class TestApplyTo:
    def test_(self):
        pass

