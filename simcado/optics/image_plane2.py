import warnings

import numpy as np
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.table import Table

from .. import utils


class ImagePlane:

    def __init__(self, header, **kwargs):

        self.meta = {"max_segment_size" : 4096*4096}
        self.meta.update(kwargs)

        if len(wcs.find_all_wcs(header=header)) == 0:
            raise ValueError("Header must have a valid WCS: {}"
                             "".format(dict(header)))

        image = np.zeros((header["NAXIS2"]+1, header["NAXIS1"]+1))
        self.hdu = fits.ImageHDU(data=image, header=header)

    def add(self, hdu_or_table, sub_pixel=False, order="bilinear"):

        if isinstance(hdu_or_table, Table):
            self.hdu.header["COMMENT"] = "Adding sources from table"
            self.hdu = add_table_to_imagehdu(hdu_or_table, self.hdu,
                                             sub_pixel=sub_pixel)
        elif isinstance(hdu_or_table, fits.ImageHDU):
            self.hdu.header["COMMENT"] = "Adding sources from table"
            self.hdu = add_imagehdu_to_imagehdu(hdu_or_table, self.hdu,
                                                order=order)

    @property
    def header(self):
        return self.hdu.header

    @property
    def data(self):
        return self.hdu.data

    @property
    def image(self):
        return self.data






def make_image_plane_header(hdu_or_table_list, pixel_scale=1*u.arcsec):
    """
    Generate a fits.Header with a WCS that covers everything in the FOV

    Parameters
    ----------
    hdu_or_table_list : list
        A list of Tables and/or ImageHDU objects for get_corner_sky_coords

    pixel_scale : astropy.Quantity
        [arcsec] The pixel scale to the projection. Default in 1 arcsec

    Returns
    -------
    header : fits.Header
        A Header containing a WCS and NAXISn values to build an ImageHDU

    See Also
    --------
    get_corner_sky_coords

    """

    pixel_scale = utils.quantify(pixel_scale, u.arcsec).to(u.deg)
    unit = pixel_scale.unit

    (x0, x1), (y0, y1) = get_sky_coords_boundaries(hdu_or_table_list)

    # old way of doing it
    # naxis1 = np.ceil((x1 - x0)*u.deg / pixel_scale).value
    # naxis2 = np.ceil((y1 - y0)*u.deg / pixel_scale).value

    the_wcs = wcs.WCS(naxis=2)
    the_wcs.wcs.crpix = [0., 0.]
    the_wcs.wcs.cdelt = [pixel_scale.value, pixel_scale.value]
    the_wcs.wcs.crval = [0., 0.]
    the_wcs.wcs.cunit = [unit, unit]
    the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    xpix, ypix = the_wcs.wcs_world2pix([x0, x1], [y0, y1], 1, ra_dec_order=True)
    naxis1 = np.round(max(xpix) - min(xpix)).astype(int)
    naxis2 = np.round(max(ypix) - min(ypix)).astype(int)

    # Add a row of padding to the image to catch outliers
    header = the_wcs.to_header()
    header["NAXIS1"] = naxis1 + 2
    header["NAXIS2"] = naxis2 + 2
    header["CRPIX1"] = -xpix[0] + 1
    header["CRPIX2"] = -ypix[0] + 1

    if naxis1 * naxis2 > 2**25:   # 2 * 4096**2
        warnings.warn("Header dimension are large: {}. Any image made from "
                      "this header will use more that >256 MB in memory"
                      "".format((naxis1, naxis2)))
    elif naxis1 * naxis2 > 2**28:
        raise MemoryError("Header dimensions are too large: {}. Any image made "
                          "from this header will use more that >8 GB in memory"
                          "".format((naxis1, naxis2)))

    return header


def get_sky_coords_boundaries(hdu_or_table_list, pixel_scale=1*u.arcsec):
    tables = [ht for ht in hdu_or_table_list if isinstance(ht, Table)]
    hdr_tables = combine_table_boundaries(tables, pixel_scale=1*u.arcsec)

    headers = [ht.header for ht in hdu_or_table_list
               if isinstance(ht, fits.ImageHDU)]
    hdr = combine_header_boundaries(headers + [hdr_tables])

    return hdr


def combine_header_boundaries(headers, pixel_scale=1*u.arcsec):
    x = []
    y = []
    for header in headers:
        w = wcs.WCS(header)
        x_foot, y_foot = w.calc_footprint(center=False).T
        x += list(x_foot)
        y += list(y_foot)

    pixel_scale = pixel_scale.to(u.deg).value
    hdr = header_from_sky_list(x, y, pixel_scale)

    return hdr


def combine_table_boundaries(tables, pixel_scale=1*u.arcsec):
    x = []
    y = []
    for table in tables:
        x_col = utils.quantity_from_table("x", table, u.arcsec).to(u.deg)
        y_col = utils.quantity_from_table("y", table, u.arcsec).to(u.deg)
        x += list(x_col.value)
        y += list(y_col.value)

    pixel_scale = pixel_scale.to(u.deg).value
    hdr = header_from_sky_list(x, y, pixel_scale)

    return hdr


def header_from_sky_list(x, y, pixel_scale):
    """x, y, pixel_scale : floats all in deg"""

    x = np.array(x)
    x[x > 180] -= 360
    x[x < -180] += 360
    print(x)

    the_wcs = wcs.WCS(naxis=2)
    the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    the_wcs.wcs.cunit = ["deg", "deg"]
    the_wcs.wcs.cdelt = [pixel_scale] * 2
    the_wcs.wcs.crval = [min(x), min(y)]
    the_wcs.wcs.crpix = [0.5, 0.5]

    hdr = the_wcs.to_header()
    hdr["NAXIS"] = 2
    hdr["NAXIS1"] = int((max(x) - min(x)) // pixel_scale) + 2
    hdr["NAXIS2"] = int((max(y) - min(y)) // pixel_scale) + 2

    return hdr







