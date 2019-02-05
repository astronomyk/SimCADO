import pytest
from copy import deepcopy

import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

import simcado.optics.image_plane as opt_imp


@pytest.fixture(scope="function")
def image_hdu_square():
    width = 100
    the_wcs = wcs.WCS(naxis=2)
    the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    the_wcs.wcs.cunit = ["arcsec", "arcsec"]
    the_wcs.wcs.cdelt = [1, 1]
    the_wcs.wcs.crval = [0, 0]
    the_wcs.wcs.crpix = [width // 2, width // 2]

    # theta = 24
    # ca, sa = np.cos(np.deg2rad(theta)), np.sin(np.deg2rad(theta))
    # the_wcs.wcs.pc = np.array([[ca, sa], [-sa, ca]])

    image = np.random.random(size=(width, width))
    hdu = fits.ImageHDU(data=image, header=the_wcs.to_header())

    return hdu


@pytest.fixture(scope="function")
def image_hdu_rect():
    width = 50
    height = 200
    the_wcs = wcs.WCS(naxis=2)
    the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    the_wcs.wcs.cunit = ["arcsec", "arcsec"]
    the_wcs.wcs.cdelt = [1, 1]
    the_wcs.wcs.crval = [0, 0]
    the_wcs.wcs.crpix = [width // 2, height // 2]

    image = np.random.random(size=(height, width))
    hdu = fits.ImageHDU(data=image, header=the_wcs.to_header())

    return hdu


@pytest.mark.usefixtures("image_hdu_square")
class TestGetSpatialExtentOfHeader:
    def test_returns_right_sky_coords_from_known_coords(self, image_hdu_square):
        xsky, ysky = opt_imp.get_corner_sky_coords_from_header(image_hdu_square.header)
        xsky = xsky*u.deg.to(u.arcsec)
        ysky = ysky*u.deg.to(u.arcsec)

        assert np.isclose(xsky[2] - xsky[1], image_hdu_square.header["NAXIS1"])
        assert np.isclose(ysky[1] - ysky[0], image_hdu_square.header["NAXIS2"])


@pytest.mark.usefixtures("image_hdu_square", "image_hdu_rect")
class TestGetImagePlaneExtentInPixels:
    def test_returns_extremes_of_two_headers(self, image_hdu_square,
                                             image_hdu_rect):
        hdu1, hdu2 = image_hdu_square, image_hdu_rect
        x_xtrm, y_xtrm = opt_imp.get_image_plane_extent_in_pixels([hdu1.header,
                                                                   hdu2.header])
        xsky = np.diff(x_xtrm)*u.deg.to(u.arcsec)
        ysky = np.diff(y_xtrm)*u.deg.to(u.arcsec)

        assert np.isclose(xsky, image_hdu_square.header["NAXIS1"])
        assert np.isclose(ysky, image_hdu_rect.header["NAXIS2"])


@pytest.mark.usefixtures("image_hdu_square", "image_hdu_rect")
class TestMakeImagePlaneHeader:
    def test_header_contains_future_naxis_pixel_sizes(self, image_hdu_square,
                                                      image_hdu_rect):
        hdr = opt_imp.make_image_plane_header([image_hdu_square.header,
                                               image_hdu_rect.header])
        assert hdr["NAXIS1"] == 100
        assert hdr["NAXIS2"] == 200

    @pytest.mark.parametrize("offset", -np.random.randint(200, 1001, 10))
    def test_header_contains_spread_out_regions(self, offset,
                                                image_hdu_square,
                                                image_hdu_rect):
        image_hdu_rect.header["CRVAL1"] += offset*u.arcsec.to(u.deg)
        hdr = opt_imp.make_image_plane_header([image_hdu_square.header,
                                               image_hdu_rect.header])
        image_width = image_hdu_square.header["NAXIS1"] // 2 + \
                      image_hdu_rect.header["NAXIS1"] // 2 + abs(offset)

        assert hdr["NAXIS1"] == image_width


class TestGetCornerSkyCoordsFromTable:
    def test_table_with_column_units_returns_right_values(self):
        x, y = [1] * u.arcmin, [1] * u.arcmin
        tbl = Table(names=["x", "y"], data=[x, y])
        xsky, ysky = opt_imp.get_corner_sky_coords_from_table(tbl)

        assert xsky[0]*u.deg == x[0].to(u.deg)
        assert ysky[0]*u.deg == y[0].to(u.deg)

    def test_table_with_meta_units_returns_right_values(self):
        x, y = [1], [1]
        tbl = Table(names=["x", "y"], data=[x, y])
        tbl.meta.update({"x_unit": u.arcmin, "y_unit": u.arcmin})
        xsky, ysky = opt_imp.get_corner_sky_coords_from_table(tbl)

        assert xsky[0] == x[0]*u.arcmin.to(u.deg)
        assert ysky[0] == y[0]*u.arcmin.to(u.deg)

    def test_table_with_default_units_returns_right_values(self):
        x, y = [60], [60]      # because default unit is arcsec
        tbl = Table(names=["x", "y"], data=[x, y])
        xsky, ysky = opt_imp.get_corner_sky_coords_from_table(tbl)

        assert pytest.approx(xsky[0] == x[0] * u.arcmin.to(u.deg))
        assert pytest.approx(ysky[0] == y[0] * u.arcmin.to(u.deg))


@pytest.mark.usefixtures("image_hdu_square")
class TestGetCornerSkyCoords:
    def test_returns_coords_for_combination_of_table_and_header(self,
                                                            image_hdu_square):
        x, y = [-100, 70] * u.arcsec, [0, 0] * u.arcsec
        tbl = Table(names=["x", "y"], data=[x, y])
        xsky, ysky = opt_imp.get_corner_sky_coords([tbl, image_hdu_square])

        assert np.all((xsky == x.to(u.deg).value))
