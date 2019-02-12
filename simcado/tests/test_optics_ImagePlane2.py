import pytest
from pytest import approx
from copy import deepcopy

import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

import simcado.optics.image_plane2 as opt_imp
import simcado.optics.image_plane_utils as impl_utils

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


PLOTS = True


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

    image = np.zeros((width, width))
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


@pytest.fixture(scope="function")
def input_table():
    x = [-10, -10, 0, 10, 10] * u.arcsec
    y = [-10, 10, 0, -10, 10] * u.arcsec
    tbl = Table(names=["x", "y"], data=[x, y])

    return tbl


@pytest.mark.usefixtures("image_hdu_square", "image_hdu_rect", "input_table")
class TestCombineTableBoundaries:
    def test_all_three_tables_are_inside_header_wcs(self, input_table):
        tbl1 = deepcopy(input_table)
        tbl2 = deepcopy(input_table)
        tbl3 = deepcopy(input_table)

        tbl2["x"] -= 25
        tbl3["y"] -= 60

        hdr = opt_imp.combine_table_boundaries([tbl1, tbl2, tbl3])

        w = wcs.WCS(hdr)
        for tbl in [tbl1, tbl2, tbl3]:
            x, y = w.wcs_world2pix(tbl["x"] / 3600., tbl["y"] / 3600., 1)
            for xi, yi in zip(x,y):
                assert xi >= 0 and xi < hdr["NAXIS1"]
                assert yi >= 0 and yi < hdr["NAXIS2"]


        if PLOTS:
            x, y = w.calc_footprint(center=False).T
            x, y = w.wcs_world2pix(x, y, 1)
            x0, y0 = w.wcs_world2pix(0, 0, 1)

            plt.plot(x, y, "b")
            plt.plot(x0, y0, "ro")
            for tbl in [tbl1, tbl2, tbl3]:
                x, y = w.wcs_world2pix(tbl["x"] / 3600., tbl["y"] / 3600., 1)
                plt.plot(x, y, "k.")

            plt.show()


@pytest.mark.usefixtures("image_hdu_square", "image_hdu_rect", "input_table")
class TestCombineImageHDUBoundaries:
    def test_all_two_imagehdus_are_inside_header_wcs(self, image_hdu_square,
                                                     image_hdu_rect):

        image_hdu_rect.header["CRVAL1"] -= 100 * u.arcsec.to(u.deg)
        image_hdu_square.header["CRVAL2"] += 100 * u.arcsec.to(u.deg)

        hdr = opt_imp.combine_header_boundaries([image_hdu_square,
                                                 image_hdu_rect])
        w = wcs.WCS(hdr)
        for im in [image_hdu_square, image_hdu_rect]:
            im_wcs = wcs.WCS(im)
            x, y = im_wcs.calc_footprint().T
            x, y = w.wcs_world2pix(x, y, 1)
            for xi, yi in zip(x, y):
                assert xi >= 0 and xi < hdr["NAXIS1"]
                assert yi >= 0 and yi < hdr["NAXIS2"]


        if PLOTS:
            for im in [image_hdu_square, image_hdu_rect]:
                im_wcs = wcs.WCS(im)
                x, y = im_wcs.calc_footprint().T
                x, y = w.wcs_world2pix(x, y, 1)
                plt.plot(x, y, "r-")

            x, y = w.calc_footprint(center=False).T
            x, y = w.wcs_world2pix(x, y, 1)
            x0, y0 = w.wcs_world2pix(0, 0, 1)

            plt.plot(x, y, "b")
            plt.plot(x0, y0, "ro")
            plt.show()


@pytest.mark.usefixtures("image_hdu_square", "image_hdu_rect", "input_table")
class TestGetSkyCoordsBoundaries:
    def test_all_5_objects_are_inside_header_wcs(self, image_hdu_square,
                                                 image_hdu_rect, input_table):

        tbl1 = deepcopy(input_table)
        tbl2 = deepcopy(input_table)
        tbl3 = deepcopy(input_table)

        tbl2["x"] -= 150
        tbl3["y"] -= 100

        image_hdu_rect.header["CRVAL1"] += 100 * u.arcsec.to(u.deg)
        image_hdu_square.header["CRVAL1"] += 0 * u.arcsec.to(u.deg)
        image_hdu_square.header["CRVAL2"] += 100 * u.arcsec.to(u.deg)

        hdr = opt_imp.get_sky_coords_boundaries([image_hdu_square, tbl1,
                                                 image_hdu_rect, tbl2, tbl3])

        w = wcs.WCS(hdr)
        for im in [image_hdu_square, image_hdu_rect]:
            im_wcs = wcs.WCS(im)
            x, y = im_wcs.calc_footprint().T
            x, y = w.wcs_world2pix(x, y, 1)
            for xi, yi in zip(x, y):
                assert xi >= 0 and xi < hdr["NAXIS1"]
                assert yi >= 0 and yi < hdr["NAXIS2"]

        for tbl in [tbl1, tbl2, tbl3]:
            x, y = w.wcs_world2pix(tbl["x"] / 3600., tbl["y"] / 3600., 1)
            for xi, yi in zip(x,y):
                assert xi >= 0 and xi < hdr["NAXIS1"]
                assert yi >= 0 and yi < hdr["NAXIS2"]


        if PLOTS:

            x, y = w.calc_footprint(center=False).T
            x, y = w.wcs_world2pix(x, y, 1)
            x0, y0 = w.wcs_world2pix(0, 0, 1)
            plt.plot(x, y, "b")
            plt.plot(x0, y0, "ro")

            for tbl in [tbl1, tbl2, tbl3]:
                x, y = w.wcs_world2pix(tbl["x"] / 3600., tbl["y"] / 3600., 1)
                plt.plot(x, y, "k.")

            for im in [image_hdu_square, image_hdu_rect]:
                im_wcs = wcs.WCS(im)
                x, y = im_wcs.calc_footprint().T
                x, y = w.wcs_world2pix(x, y, 1)
                plt.plot(x, y, "r-")

            plt.show()
