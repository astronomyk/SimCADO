import pytest
from pytest import approx

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


def _table_source():
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


def _image_source(dx=0, dy=0, angle=0, weight=1):
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

    im = np.random.random(size=(n+1, n+1)) * 1E-10
    im[0, n] += 5
    im[n, 0] += 5
    im[n//2, n//2] += 10

    im_hdu = fits.ImageHDU(data=im, header=im_wcs.to_header())
    im_hdu.header["SPEC_REF"] = 0
    im_source = Source(image_hdu=im_hdu, spectra=specs)

    angle = angle * np.pi / 180
    im_source.fields[0].header["CRVAL1"] += dx * u.arcsec.to(u.deg)
    im_source.fields[0].header["CRVAL2"] += dy * u.arcsec.to(u.deg)
    im_source.fields[0].header["PC1_1"] = np.cos(angle)
    im_source.fields[0].header["PC1_2"] = np.sin(angle)
    im_source.fields[0].header["PC2_1"] = -np.sin(angle)
    im_source.fields[0].header["PC2_2"] = np.cos(angle)
    im_source.fields[0].data *= weight

    return im_source


def _combined_source(im_angle=0, dx=[0, 0, 0], dy=[0, 0, 0], weight=[1, 1, 1]):
    tblsrc1 = _table_source()

    tblsrc2 = _table_source()
    tblsrc2.fields[0]["x"] += dx[0]
    tblsrc2.fields[0]["y"] += dy[0]
    tblsrc2.fields[0]["weight"] *= weight[0]

    tblsrc3 = _table_source()
    tblsrc3.fields[0]["x"] += dx[1]
    tblsrc3.fields[0]["y"] += dy[1]
    tblsrc3.fields[0]["weight"] *= weight[1]

    imsrc = _image_source(dx[2], dy[2], im_angle, weight[2])

    src = tblsrc1 + tblsrc2 + tblsrc3 + imsrc

    return src


@pytest.fixture(scope="function")
def table_source():
    return _table_source()


@pytest.fixture(scope="function")
def image_source():
    return _image_source()


@pytest.fixture(scope="function")
def combined_source():
    return _combined_source()


@pytest.fixture(scope="function")
def basic_fov_header():
    w, h = 150, 150
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

    def test_throws_error_if_no_wcs_in_header(self):
        with pytest.raises(ValueError):
            fov.FieldOfView(fits.Header(), (1, 2) * u.um)


@pytest.mark.usefixtures("basic_fov_header", "table_source", "image_source")
class TestFieldOfViewExtractFrom:
    def test_creates_single_combined_from_multiple_tables(self, table_source,
                                                          basic_fov_header):
        src = table_source + table_source
        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um)
        the_fov.extract_from(src)
        assert len(the_fov.fields) == 1
        assert isinstance(the_fov.fields[0], Table)

    def test_creates_single_combined_from_multiple_images(self, image_source,
                                                          basic_fov_header):
        src = image_source + image_source
        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um)
        the_fov.extract_from(src)
        assert len(the_fov.fields) == 1
        assert isinstance(the_fov.fields[0], fits.ImageHDU)

    def test_creates_two_fields_for_tables_and_images(self, basic_fov_header,
                                                      image_source,
                                                      table_source):
        src = image_source + table_source
        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um)
        the_fov.extract_from(src)
        assert len(the_fov.fields) == 2
        assert isinstance(the_fov.fields[0], Table)
        assert isinstance(the_fov.fields[1], fits.ImageHDU)

    def test_ignores_fields_outside_fov_boundary(self, basic_fov_header):

        src = _combined_source(dx=[200, 200, 200])
        src.fields[0]["x"] += 200

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um)
        the_fov.extract_from(src)

        assert len(the_fov.fields) == 0


@pytest.mark.usefixtures("basic_fov_header")
class TestFieldOfViewView:
    def test_views_with_only_image(self, basic_fov_header):
        src = _image_source()
        flux = src.photons_in_range(1*u.um, 2*u.um).value
        orig_sum = np.sum(src.fields[0].data * flux)

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um)
        the_fov.extract_from(src)
        view = the_fov.view()

        assert np.sum(view) == approx(orig_sum)

        if PLOTS:
            plt.imshow(src.fields[0].data.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

            plt.imshow(view.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

    def test_views_with_rotated_image(self, basic_fov_header):
        src = _image_source(angle=30)
        flux = src.photons_in_range(1*u.um, 2*u.um).value
        orig_sum = np.sum(src.fields[0].data) * flux

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um)
        the_fov.extract_from(src)
        view = the_fov.view()

        assert np.sum(view) == approx(orig_sum)

        if PLOTS:
            plt.imshow(src.fields[0].data.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

            plt.imshow(view.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

    def test_views_with_only_table(self, basic_fov_header):
        src = _table_source()
        fluxes = src.photons_in_range(1*u.um, 2*u.um)
        phs = fluxes[src.fields[0]["ref"]] * src.fields[0]["weight"]

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2) * u.um)
        the_fov.extract_from(src)
        view = the_fov.view()

        assert np.sum(view) == np.sum(phs).value - 4

        if PLOTS:
            plt.imshow(view.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()

    def test_views_with_tables_and_images(self, basic_fov_header):
        src = _combined_source(im_angle=0, weight=[1, 1, 1],
                               dx=[-2, 0, 0], dy=[2, 2, 0])
        ii = src.fields[3].header["SPEC_REF"]
        flux = src.photons_in_range(1*u.um, 2*u.um, indexes=[ii]).value
        orig_sum = np.sum(src.fields[3].data) * flux

        the_fov = fov.FieldOfView(basic_fov_header, (1, 2)*u.um)
        the_fov.extract_from(src)
        view = the_fov.view()

        assert np.sum(the_fov.fields[0]["flux"]) == approx(36)
        assert np.sum(the_fov.fields[1].data) == approx(orig_sum)
        assert np.sum(view) == approx(orig_sum + 36)

        if PLOTS:
            plt.imshow(view.T, origin="lower", norm=LogNorm())
            plt.colorbar()
            plt.show()


@pytest.mark.usefixtures("basic_fov_header", "image_source", "table_source")
class TestIsFieldInFOV:
    def test_returns_true_for_table_inside(self, basic_fov_header,
                                           table_source):
        assert fov.is_field_in_fov(basic_fov_header, table_source.fields[0])

    def test_returns_false_for_table_outside(self, basic_fov_header,
                                             table_source):
        table_source.fields[0]["x"] += 200
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
        image_source.fields[0].header["CRVAL1"] += 200*u.arcsec.to(u.deg)
        assert not fov.is_field_in_fov(basic_fov_header, image_source.fields[0])

    def test_returns_true_for_image_in_corner(self, basic_fov_header,
                                              image_source):
        image_source.fields[0].header["CRVAL1"] += 10*u.arcsec.to(u.deg)
        image_source.fields[0].header["CRVAL2"] -= 10*u.arcsec.to(u.deg)
        assert fov.is_field_in_fov(basic_fov_header, image_source.fields[0])


class TestMakeFluxTable:
    def test_flux_in_equals_flux_out(self):
        pass


class TestCombineTableFields:
    def test_flux_in_equals_flux_out(self):
        pass


class TestCombineImageHDUFields:
    def test_flux_in_equals_flux_out(self):
        pass


class TestHasWcsKeys:
    def test_fails_if_header_does_not_have_all_keys(self):
        assert not fov.has_needed_keywords(fits.Header())

    def test_passes_if_header_does_have_all_keys(self):
        hdr = wcs.WCS().to_header()
        hdr["NAXIS1"] = 100
        assert fov.has_needed_keywords(hdr)
