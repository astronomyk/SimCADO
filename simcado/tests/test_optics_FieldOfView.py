import pytest
from pytest import approx
from copy import deepcopy

import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

from synphot import SourceSpectrum
from synphot.models import Empirical1D

from simcado.optics import fov
from simcado.source.source2 import Source

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

PLOTS = False


@pytest.fixture(scope="function")
def table_source():
    n = 100
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=4 * np.ones(n) * unit),
             SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n) * unit),
             SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n)[::-1] * unit)]
    tbl = Table(names=["x", "y", "ref", "weight"],
                data=[[5,  0, -5,  0]*u.arcsec,
                      [5, -10, 5,  0] * u.arcsec,
                      [2,  0,  1,  0],
                      [1,  1,  1,  2]])
    tbl_source = Source(table=tbl, spectra=specs)

    return tbl_source


@pytest.fixture(scope="function")
def image_source():
    n = 50
    unit = u.Unit("ph s-1 m-2 um-1")
    wave = np.linspace(0.5, 2.5, n) * u.um
    specs = [SourceSpectrum(Empirical1D, points=wave,
                            lookup_table=np.linspace(0, 4, n) * unit)]

    n = 50
    im_wcs = wcs.WCS(naxis=2)
    im_wcs.wcs.cunit = [u.arcsec, u.arcsec]
    im_wcs.wcs.cdelt = [0.2, 0.2]
    im_wcs.wcs.crval = [0, 0]
    im_wcs.wcs.crpix = [n//2, n//2]
    im_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    im = np.ones((n+1, n+1)) * 1E-11
    im[0, n] += 5
    im[n, 0] += 5
    im[n//2, n//2] += 10

    im_hdu = fits.ImageHDU(data=im, header=im_wcs.to_header())
    im_hdu.header["SPEC_REF"] = 0
    im_source = Source(image_hdu=im_hdu, spectra=specs)

    return im_source


@pytest.fixture(scope="function")
def basic_fov_header():
    w, h = 100, 100
    fovwcs = wcs.WCS(naxis=2)
    fovwcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    fovwcs.wcs.cdelt = [0.1, 0.1]
    fovwcs.wcs.cunit = ["arcsec", "arcsec"]
    fovwcs.wcs.crval = [0, 0]
    fovwcs.wcs.crpix = [w / 2, h / 2]

    fovhdr = fovwcs.to_header()
    fovhdr["NAXIS"] = 2
    fovhdr["NAXIS1"] = w
    fovhdr["NAXIS2"] = h

    return fovhdr


@pytest.mark.usefixtures("basic_fov_header")
class TestFieldOfViewInit:
    def test_initialises_with_nothing_raise_error(self):
        with pytest.raises(TypeError):
            fov.FieldOfView()

    def test_initialises_with_header_and_waverange(self, basic_fov_header):
        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um)
        assert isinstance(the_fov, fov.FieldOfView)


@pytest.mark.usefixtures("basic_fov_header", "image_source", "table_source")
class TestFieldOfViewExtractFrom:
    def test_temp_extract_from(self, basic_fov_header, table_source,
                               image_source):
        tblsrc1 = deepcopy(table_source)
        tblsrc2 = deepcopy(table_source)
        tblsrc2.fields[0]["x"] += 8   # puts it outside
        tblsrc3 = deepcopy(table_source)
        tblsrc3.fields[0]["y"] += 12   # still inside
        tblsrc3.fields[0]["weight"] *= 1

        image_source.fields[0].header["CRVAL1"] += 1*u.arcsec.to(u.deg)

        src = tblsrc1 + image_source + tblsrc2 + tblsrc3
        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um)

        the_fov.extract_from(src)

        if PLOTS:
            plt.imshow(the_fov.fields[1].data.T, origin="lower")
            plt.colorbar()
            plt.show()









@pytest.mark.usefixtures("basic_fov_header", "image_source", "table_source")
class TestIsFieldInFOV:
    def test_returns_true_for_table_inside(self, basic_fov_header,
                                           table_source):
        assert fov.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_false_for_table_outside(self, basic_fov_header,
                                             table_source):
        table_source.fields[0]["x"] += 12
        assert not fov.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_true_for_table_corner_inside(self, basic_fov_header,
                                                  table_source):
        table_source.fields[0]["x"] += 9
        table_source.fields[0]["y"] += 14
        assert fov.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_true_for_image_inside(self, basic_fov_header,
                                           image_source):
        assert fov.is_field_in_fov(basic_fov_header, image_source.fields[0])

    def test_returns_false_for_image_outside(self, basic_fov_header,
                                             image_source):
        image_source.fields[0].header["CRVAL1"] += 11*u.arcsec.to(u.deg)
        assert not fov.is_field_in_fov(basic_fov_header, image_source.fields[0])

    def test_returns_true_for_image_in_corner(self, basic_fov_header,
                                              image_source):
        image_source.fields[0].header["CRVAL1"] += 10*u.arcsec.to(u.deg)
        image_source.fields[0].header["CRVAL2"] -= 10*u.arcsec.to(u.deg)
        assert fov.is_field_in_fov(basic_fov_header, image_source.fields[0])





