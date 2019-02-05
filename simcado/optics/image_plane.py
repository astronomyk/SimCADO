import warnings

import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.table import Table

from .. import utils


def get_corner_sky_coords(positions):

    if isinstance(positions, (list, tuple)):
        xlist, ylist = [], []
        for position in positions:
            x, y = get_corner_sky_coords(position)
            xlist += list(x)
            ylist += list(y)

        xsky = [np.min(xlist), np.max(xlist)]
        ysky = [np.min(ylist), np.max(ylist)]

    else:
        if isinstance(positions, Table):
            xsky, ysky = get_corner_sky_coords_from_table(positions)
        elif isinstance(positions, fits.ImageHDU):
            xsky, ysky = get_corner_sky_coords_from_header(positions.header)
        else:
            xsky, ysky = [], []

    return xsky, ysky


def get_corner_sky_coords_from_table(table):
    if not isinstance(table, Table):
        raise ValueError("table must be astropy.Table: {}".format(type(table)))

    x = utils.quantity_from_table("x", table, default_unit=u.arcsec).to(u.deg)
    y = utils.quantity_from_table("y", table, default_unit=u.arcsec).to(u.deg)
    x[x > 180. * u.deg] -= 360. * u.deg

    xsky = [np.min(x).value, np.max(x).value]
    ysky = [np.min(y).value, np.max(y).value]

    return xsky, ysky


def get_corner_sky_coords_from_header(header):
    if not isinstance(header, fits.Header):
        raise ValueError("header must be astropy.Header: {}"
                         "".format(type(header)))
    naxis1, naxis2 = header["NAXIS1"], header["NAXIS2"]
    xpix, ypix = [0, 0, naxis1, naxis1], [0, naxis2, naxis2, 0]

    hdr_wcs = wcs.WCS(header)
    xsky, ysky = hdr_wcs.wcs_pix2world(xpix, ypix, 1)
    xsky[xsky > 180.] -= 360.

    return xsky, ysky


def get_image_plane_extent_in_pixels(headers):

    xsky, ysky = [], []
    for header in headers:
        naxis1, naxis2 = header["NAXIS1"], header["NAXIS2"]
        xpix, ypix = [0, 0, naxis1, naxis1], [0, naxis2, naxis2, 0]

        hdr_wcs = wcs.WCS(header)
        x, y = hdr_wcs.wcs_pix2world(xpix, ypix, 1)
        x[x > 180.] -= 360.
        xsky += list(x)
        ysky += list(y)

    return (np.min(xsky), np.max(xsky)), \
           (np.min(ysky), np.max(ysky))


def make_image_plane_header(headers, pixel_scale=1*u.arcsec):

    pixel_scale = pixel_scale.to(u.deg)
    unit = pixel_scale.unit

    (x0, x1), (y0, y1) = get_image_plane_extent_in_pixels(headers)
    # naxis1 = np.ceil((x1 - x0)*u.deg / pixel_scale).value
    # naxis2 = np.ceil((y1 - y0)*u.deg / pixel_scale).value

    the_wcs = wcs.WCS(naxis=2)
    the_wcs.wcs.crpix = [0., 0.]
    the_wcs.wcs.cdelt = [pixel_scale.value, pixel_scale.value]
    the_wcs.wcs.crval = [0., 0.]
    the_wcs.wcs.cunit = [unit, unit]
    the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    xpix, ypix = the_wcs.wcs_world2pix([x0, x1], [y0, y1], 1)
    naxis1 = np.round(xpix[1] - xpix[0]).astype(int)
    naxis2 = np.round(ypix[1] - ypix[0]).astype(int)

    hdr = the_wcs.to_header()
    hdr["NAXIS1"] = naxis1
    hdr["NAXIS2"] = naxis2
    hdr["CRPIX1"] = -xpix[0]
    hdr["CRPIX2"] = -ypix[0]

    if naxis1 * naxis2 > 2**25:   # 2 * 4096**2
        warnings.warn("Header dimension are large: {}. Any image made from"
                      "this header will use more that >256 MB in memory"
                      "".format((naxis1, naxis2)))

    return hdr


class ImagePlane:
    def __init__(self, axis_dims=None, pixel_scale=None, **kwargs):

        self.meta = {"max_segment_size"}
        self.meta.update(kwargs)




