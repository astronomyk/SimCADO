import warnings

import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.table import Table

from reproject import reproject_interp

from .. import utils


class ImagePlane:
    """
    A class to act as a canvas onto which to project `Source` images or tables

    Parameters
    ----------
    header : `fits.Header`
        Must contain a valid WCS

    kwargs
    ------
    max_segment_size : int
        [pixel] Default is 4096**2. Maximum amount of pixels per segment.
        Used in optimising memory usage and computation overheads

    .. todo: Write the code to deal with a canvas larger than max_segment_size

    Examples
    --------
    ::

        from astropy.table import Table
        from simcado.optics import image_plane as impl

        my_point_source_table = Table(names=["x", "y", "flux"],
                                      data=[(0,  1,  2)*u.arcsec,
                                            (0, -5, 10)*u.arcsec,
                                            (100,50,25)*u.ph/u.s])

        hdr = impl.make_image_plane_header([my_point_source_table],
                                           pixel_scale=0.1*u.arcsec)
        img_plane = impl.ImagePlane(hdr)
        img_plane.add(my_point_source_table)

        print(img_plane.image)

    """

    def __init__(self, header, **kwargs):

        self.meta = {"max_segment_size" : 4096*4096}
        self.meta.update(kwargs)

        if len(wcs.find_all_wcs(header=header)) == 0:
            raise ValueError("Header must have a valid WCS: {}"
                             "".format(dict(header)))

        image = np.zeros((header["NAXIS1"]+1, header["NAXIS2"]+1))
        self.hdu = fits.ImageHDU(data=image, header=header)

    def add(self, hdu_or_table, sub_pixel=False):
        """
        Add a projection of an image or table sources to the canvas

        .. note::
          If a Table is provides, it must include the following columns:
          `x`, `y`, and `flux`.

          Units for the columns should be provided in the
          <Table>.unit attribute or as an entry in the table's meta dictionary
          using this syntax: <Table>.meta["<colname>_unit"] = <unit>.

          For example::

              tbl["x"].unit = u.arcsec
              tbl.meta["x_unit"] = "deg"

          If no units are given, default units will be assumed. These are:

          - `x`, `y`: `arcsec`
          - `flux` : `ph / s / pix`

        Parameters
        ----------
        hdu_or_table : `fits.ImageHDU` or `astropy.Table`
            The input to be projected onto the image plane. See above.

        sub_pixel : bool, optional
            Default is False. Dictates if point sources should be projected with
            sub-pixel shifts or not. Accounting for sub-pixel shifts is approx.
            5x slower.

        """

        if isinstance(hdu_or_table, Table):
            self.hdu.header["COMMENT"] = "Adding sources from table"
            self.hdu = add_table_to_imagehdu(hdu_or_table, self.hdu,
                                             sub_pixel=sub_pixel)
        elif isinstance(hdu_or_table, fits.ImageHDU):
            self.hdu.header["COMMENT"] = "Adding sources from table"
            self.hdu = add_imagehdu_to_imagehdu(hdu_or_table, self.hdu)

    @property
    def data(self):
        return self.hdu.data

    @property
    def header(self):
        return self.hdu.header

    @property
    def image(self):
        return self.data


def get_corner_sky_coords(hdus_tables):
    """
    Get the sky coordinates for the corners of a position object

    Parameters
    ----------
    hdus_tables : list, fits.ImageHDU, astropy.Table

    Returns
    -------
    x_edges, y_edges : list
        [x_min, x_max], [y_min, y_max]

    """

    if isinstance(hdus_tables, (list, tuple)):
        xlist, ylist = [], []
        for hdu_or_table in hdus_tables:
            x, y = get_corner_sky_coords(hdu_or_table)
            xlist += list(x)
            ylist += list(y)

        x_edges = [np.min(xlist), np.max(xlist)]
        y_edges = [np.min(ylist), np.max(ylist)]

    else:
        if isinstance(hdus_tables, Table):
            x_edges, y_edges = get_corner_sky_coords_from_table(hdus_tables)
        elif isinstance(hdus_tables, fits.ImageHDU):
            hdr = hdus_tables.header
            x_edges, y_edges = get_corner_sky_coords_from_header(hdr)
        else:
            x_edges, y_edges = [], []

    return x_edges, y_edges


def get_corner_sky_coords_from_table(table):
    """
    Returns the extreme x, y values found in an astropy.Table

    Parameters
    ----------
    table : astropy.Table
        Must contain columns "x", "y"  preferably with units in the .unit
        attribute of each column. Alternatively include "x_unit" and "y_unit"
        in the <table>.meta dictionary

    Returns
    -------
    x_edges, y_edges : list
        [x_min, x_max], [y_min, y_max]

    """

    if not isinstance(table, Table):
        raise ValueError("table must be astropy.Table: {}".format(type(table)))

    x = utils.quantity_from_table("x", table, default_unit=u.arcsec).to(u.deg)
    y = utils.quantity_from_table("y", table, default_unit=u.arcsec).to(u.deg)
    x[x > 180. * u.deg] -= 360. * u.deg

    x_edges = [np.min(x).value, np.max(x).value]
    y_edges = [np.min(y).value, np.max(y).value]

    return x_edges, y_edges


def get_corner_sky_coords_from_header(header):
    """
    Returns the extreme x, y values found in an fits.Header object

    Parameters
    ----------
    header : fits.Header
        Must contain a valid WCS, i.e. CTYPEn, CDELTn, CRPIXn, CRVALn, CUNITn

    Returns
    -------
    x_edges, y_edges : list
        [x_min, x_max], [y_min, y_max]

    """

    if not isinstance(header, fits.Header):
        raise ValueError("header must be astropy.Header: {}"
                         "".format(type(header)))

    naxis1, naxis2 = header["NAXIS1"], header["NAXIS2"]
    xpix, ypix = [0, 0, naxis1, naxis1], [0, naxis2, naxis2, 0]

    hdr_wcs = wcs.WCS(header)
    x_edges, y_edges = hdr_wcs.wcs_pix2world(xpix, ypix, 1)
    x_edges[x_edges > 180.] -= 360.

    return x_edges, y_edges


def get_image_plane_extent_in_pixels(headers):
    """ Depreciated by get_corner_sky_coords()"""

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


def make_image_plane_header(hdu_or_table_list, pixel_scale=1 * u.arcsec):
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

    (x0, x1), (y0, y1) = get_corner_sky_coords(hdu_or_table_list)
    # old way of doing it
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

    header = the_wcs.to_header()
    header["NAXIS1"] = naxis1
    header["NAXIS2"] = naxis2
    header["CRPIX1"] = -xpix[0]
    header["CRPIX2"] = -ypix[0]

    if naxis1 * naxis2 > 2**25:   # 2 * 4096**2
        warnings.warn("Header dimension are large: {}. Any image made from "
                      "this header will use more that >256 MB in memory"
                      "".format((naxis1, naxis2)))
    elif naxis1 * naxis2 > 2**28:
        raise MemoryError("Header dimensions are too large: {}. Any image made "
                          "from this header will use more that >8 GB in memory"
                          "".format((naxis1, naxis2)))

    return header


def add_table_to_imagehdu(table, canvas_hdu, sub_pixel=True):
    """
    Add sources from an astropy.Table to the image of an fits.ImageHDU

    Parameters
    ----------
    table : astropy.Table
        Must contain the columns "x", "y", "flux" with the units in the column
        attribute .unit, or in the table.meta dictionary as "<colname>_unit"

    canvas_hdu : fits.ImageHDU
        The ImageHDU onto which the table sources should be projected.
        This must include a valid WCS

    sub_pixel : bool, optional
        Default is True. If True, sub-pixel shifts of sources will be taken into
        account when projecting onto the canvas pixel grid. This takes about 5x
        longer than ignoring the sub-pixel shifts

    Returns
    -------
    canvas_hdu : fits.ImageHDU

    """

    if len(wcs.find_all_wcs(canvas_hdu.header)) == 0:
        raise ValueError("canvas_hdu must include a WCS")

    x = utils.quantity_from_table("x", table, default_unit=u.arcsec).to(u.deg)
    y = utils.quantity_from_table("y", table, default_unit=u.arcsec).to(u.deg)
    f = utils.quantity_from_table("flux", table, default_unit=u.Unit("ph s-1"))
    x[x > 180. * u.deg] -= 360. * u.deg

    canvas_wcs = wcs.WCS(canvas_hdu)
    canvas_hdu.header["comment"] = "x_unit : {}".format(x.unit)
    canvas_hdu.header["comment"] = "y_unit : {}".format(y.unit)
    canvas_hdu.header["comment"] = "flux_unit : {}".format(f.unit)

    xpix, ypix = canvas_wcs.wcs_world2pix(x, y, 1)

    if sub_pixel is True:
        canvas_hdu.header["comment"] = "Adding {} sub-pixel sources" \
                                       "".format(len(f))
        for ii in range(len(xpix)):
            xx, yy, fracs = sub_pixel_fractions(xpix[ii], ypix[ii])
            for x, y, frac in zip(xx, yy, fracs):
                canvas_hdu.data[x, y] += frac * f[ii].value
    else:
        canvas_hdu.header["comment"] = "Adding {} int-pixel sources" \
                                       "".format(len(f))
        xpix = xpix.astype(int)
        ypix = ypix.astype(int)
        for ii in range(len(xpix)):
            canvas_hdu.data[xpix[ii], ypix[ii]] += f[ii].value

    return canvas_hdu


def add_imagehdu_to_imagehdu(image_hdu, canvas_hdu, order="bilinear"):
    """
    Re-project one ``fits.ImageHDU`` onto another ``fits.ImageHDU``

    Parameters
    ----------
    image_hdu : fits.ImageHDU
        The ``ImageHDU`` which will be reprojected onto `canvas_hdu`


    canvas_hdu : fits.ImageHDU
        The ``ImageHDU`` onto which the table sources should be projected.
        This must include a valid WCS

    order : str, optional
        Default is ``bilinear``. The resampling scheme used by
        ``reproject_interp``. See ``reproject.reproject_interp()`` for the list
        of options

    Returns
    -------
    canvas_hdu : fits.ImageHDU

    """

    new_im, mask = reproject_interp(image_hdu, canvas_hdu.header, order=order)
    new_im[mask > 0] *= np.sum(image_hdu.data) / np.sum(new_im[mask > 0])
    canvas_hdu.data[mask > 0] += new_im[mask > 0]

    return canvas_hdu


def sub_pixel_fractions(x, y):
    """
    Makes a list of pixel coordinates and weights to reflect sub-pixel shifts

    A point source which isn't centred on a pixel can be modelled by a centred
    PSF convolved with a shifted delta function. A fraction of the delta
    function in moved into each of the adjoining pixels. For example, a star
    at ``(x,y)=(0.2, 0.2)`` would be represented by a following pixel weights::

        ---------------
        | 0.16 | 0.04 |
        ---------------
        | 0.64 | 0.16 |
        ---------------

    where (0,0) is the centre of the bottom-left pixel

    Given (x,y) pixel coordinates, this function returns the fractions of flux
    that should go into the surrounding pixels, as well as the coordinates of
    those neighbouring pixels.

    Parameters
    ----------
    x, y : float

    Returns
    -------
    x_pix, y_pix, fracs : list of (int, int, float)
        The x and y pixel coordinates and their corresponding flux fraction

    """

    x0, dx = divmod(x, 1)
    y0, dy = divmod(y, 1)

    xi0 = int(x0)
    xi1 = xi0 + bool(dx)
    yi0 = int(y0)
    yi1 = yi0 + bool(dy)

    f00 = (1. - dx) * (1. - dy)
    f01 = (1. - dx) * dy
    f10 = dx * (1 - dy)
    f11 = dx * dy

    x_pix = [xi0, xi1, xi0, xi1]
    y_pix = [yi0, yi0, yi1, yi1]
    fracs = [f00, f10, f01, f11]

    return x_pix, y_pix, fracs
