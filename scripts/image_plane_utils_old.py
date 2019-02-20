import warnings

import numpy as np   # must be above 1.13 for the np.num_to_nan(copy=...) to work

from astropy import units as u, wcs
from astropy.io import fits
from astropy.table import Table
from reproject import reproject_interp

from .. import utils


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

    canvas_hdu.header["comment"] = "x_unit : {}".format(x.unit)
    canvas_hdu.header["comment"] = "y_unit : {}".format(y.unit)
    canvas_hdu.header["comment"] = "flux_unit : {}".format(f.unit)

    canvas_wcs = wcs.WCS(canvas_hdu)
    xpix, ypix = canvas_wcs.wcs_world2pix(x, y, 1)

    naxis1 = canvas_hdu.header["NAXIS1"]
    naxis2 = canvas_hdu.header["NAXIS2"]
    # Weird FITS / astropy behaviour. Axis1 == y, Axis2 == x.
    # Also occasionally 0 is returned as ~ -1e-11
    eps = -1e-7
    mask = (xpix >= eps) * (xpix < naxis2) * (ypix >= eps) * (ypix < naxis1)

    if sub_pixel is True:
        canvas_hdu.header["comment"] = "Adding {} sub-pixel sources" \
                                       "".format(len(f))
        for ii in range(len(xpix)):
            if mask[ii]:
                xx, yy, fracs = sub_pixel_fractions(xpix[ii], ypix[ii])
                for x, y, frac in zip(xx, yy, fracs):
                    canvas_hdu.data[x, y] += frac * f[ii].value
    else:
        canvas_hdu.header["comment"] = "Adding {} int-pixel sources" \
                                       "".format(len(f))
        xpix = xpix.astype(int)
        ypix = ypix.astype(int)
        for ii in range(len(xpix)):
            if mask[ii]:
                canvas_hdu.data[xpix[ii], ypix[ii]] += f[ii].value

    return canvas_hdu


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
        [deg] : [x_min, x_max], [y_min, y_max] on sky coordinates

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

    xsky, ysky = wcs.WCS(header).calc_footprint(center=False).T
    xsky[xsky > 180] -= 360
    x_edges = [min(xsky), max(xsky)]
    y_edges = [min(ysky), max(ysky)]

    return x_edges, y_edges


def make_image_plane_header(hdu_or_table_list, pixel_scale=1*u.arcsec):
    """
    Generate a fits.Header with a WCS that covers everything in the FOV

    Parameters
    ----------
    hdu_or_table_list : list
        A list of Tables and/or ImageHDU py_objects for get_corner_sky_coords

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


def overlay_image(small_im, big_im, coords, mask=None, sub_pixel=False):
    """
    Overlay small_im on top of big_im at the position specified by coords

    ``small_im`` will be centred at ``coords``

    Adapted from:
    ``https://stackoverflow.com/questions/14063070/
        overlay-a-smaller-image-on-a-larger-image-python-opencv``

    """

    # TODO - Add in a catch for sub-pixel shifts
    if sub_pixel:
        raise NotImplementedError

    x, y = np.array(coords, dtype=int) - np.array(small_im.shape) // 2

    # Image ranges
    x1, x2 = max(0, x), min(big_im.shape[0], x + small_im.shape[0])
    y1, y2 = max(0, y), min(big_im.shape[1], y + small_im.shape[1])

    # Overlay ranges
    x1o, x2o = max(0, -x), min(small_im.shape[0], big_im.shape[0] - x)
    y1o, y2o = max(0, -y), min(small_im.shape[1], big_im.shape[1] - y)

    # Exit if nothing to do
    if y1 >= y2 or x1 >= x2 or y1o >= y2o or x1o >= x2o:
        return big_im

    if mask is None:
        big_im[x1:x2, y1:y2] += small_im[x1o:x2o, y1o:y2o]
    else:
        mask = mask[x1o:x2o, y1o:y2o].astype(bool)
        big_im[x1:x2, y1:y2][mask] += small_im[x1o:x2o, y1o:y2o][mask]

    return big_im


def overlay_imagehdu_on_imagehdu(image_hdu, canvas_hdu, mask=None):
    wcs_image = wcs.WCS(image_hdu)
    wcs_canvas = wcs.WCS(canvas_hdu)

    xcen_im = image_hdu.header["NAXIS1"] // 2
    ycen_im = image_hdu.header["NAXIS2"] // 2

    ysky0, xsky0 = wcs_image.wcs_pix2world(xcen_im, ycen_im, 1)
    xpix0, ypix0 = wcs_canvas.wcs_world2pix(xsky0, ysky0, 1)
    canvas_hdu.data = overlay_image(image_hdu.data, canvas_hdu.data,
                                    coords=(xpix0, ypix0), mask=mask)

    return canvas_hdu


def make_canvas_header_for_reproject(in_imagehdu, pixel_scale=4*u.mas):
    wcs_image = wcs.WCS(in_imagehdu)
    footprint = wcs_image.calc_footprint().T
    xsky_edges, ysky_edges = footprint[0], footprint[1]

    xsky_min, ysky_min = min(xsky_edges), min(ysky_edges)

    pixel_scale = pixel_scale.to(u.deg).value
    wcs_canvas = wcs.WCS(naxis=2)
    wcs_canvas.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs_canvas.wcs.cunit = ["deg", "deg"]
    wcs_canvas.wcs.cdelt = [pixel_scale, pixel_scale]
    wcs_canvas.wcs.crval = [xsky_min, ysky_min]
    wcs_canvas.wcs.crpix = [0, 0]

    xpix, ypix = wcs_canvas.wcs_world2pix(xsky_edges, ysky_edges, 1)
    hdr_canvas = wcs_canvas.to_header()
    hdr_canvas["NAXIS"] = 2
    hdr_canvas["NAXIS1"] = int(np.ceil(max(xpix) - min(xpix)))
    hdr_canvas["NAXIS2"] = int(np.ceil(max(ypix) - min(ypix)))

    return hdr_canvas


def add_imagehdu_to_imagehdu(image_hdu, canvas_hdu, order="bilinear"):
    if isinstance(image_hdu.data, u.Quantity):
        unit = image_hdu.data.unit
        image_hdu.data = image_hdu.data.value

    image_hdu.header["CRVAL1"] *= -1
    image_hdu.header["CRVAL1"] += 180.
    canvas_hdu.header["CRVAL1"] += 180.

    pixel_scale = canvas_hdu.header["CDELT1"] * u.deg
    hdr_reproj = make_canvas_header_for_reproject(image_hdu,
                                                  pixel_scale=pixel_scale)
    new_im, mask = reproject_interp(image_hdu, hdr_reproj, order=order)
    new_im = np.nan_to_num(new_im, copy=False)
    new_im[mask > 0] *= np.sum(image_hdu.data) / np.sum(new_im[mask > 0])
    new_im_hdu = fits.ImageHDU(data=new_im, header=hdr_reproj)

    canvas_hdu = overlay_imagehdu_on_imagehdu(new_im_hdu, canvas_hdu, mask=mask)

    canvas_hdu.header["CRVAL1"] -= 180.
    image_hdu.header["CRVAL1"] -= 180.
    image_hdu.header["CRVAL1"] *= -1

    return canvas_hdu


# def add_imagehdu_to_imagehdu(image_hdu, canvas_hdu, order="bilinear"):
#     """
#     Re-project one ``fits.ImageHDU`` onto another ``fits.ImageHDU``
#
#     Parameters
#     ----------
#     image_hdu : fits.ImageHDU
#         The ``ImageHDU`` which will be reprojected onto `canvas_hdu`
#
#
#     canvas_hdu : fits.ImageHDU
#         The ``ImageHDU`` onto which the table sources should be projected.
#         This must include a valid WCS
#
#     order : str, optional
#         Default is ``bilinear``. The resampling scheme used by
#         ``reproject_interp``. See ``reproject.reproject_interp()`` for the list
#         of options
#
#     Returns
#     -------
#     canvas_hdu : fits.ImageHDU
#
#     """
#
#     if isinstance(image_hdu.data, u.Quantity):
#         unit = image_hdu.data.unit
#         image_hdu.data = image_hdu.data.value
#
#     new_im, mask = reproject_interp(image_hdu, canvas_hdu.header, order=order)
#     new_im = np.nan_to_num(new_im, copy=False)
#
#
#     new_wcs = wcs.WCS(canvas_hdu)
#     x, y = new_wcs.wcs_world2pix([0], [0], 1)
#     import matplotlib.pyplot as plt
#     from matplotlib.colors import LogNorm
#     plt.imshow(new_im.T, origin="lower", norm=LogNorm())
#     plt.scatter(x, y, c="r")
#     plt.show()
#
#     # this won't work when image_hdu is larger than canvas_hdu
#     if np.prod(canvas_hdu.data.shape) > np.prod(image_hdu.data.shape):
#         new_im[mask > 0] *= np.sum(image_hdu.data) / np.sum(new_im[mask > 0])
#
#     canvas_hdu.data[mask > 0] += new_im[mask > 0]
#
#     return canvas_hdu

