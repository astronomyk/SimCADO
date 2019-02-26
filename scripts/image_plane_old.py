import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.table import Table

from .image_plane_utils_old import add_table_to_imagehdu, add_imagehdu_to_imagehdu


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

        image = np.zeros((header["NAXIS2"]+1, header["NAXIS1"]+1))
        self.hdu = fits.ImageHDU(data=image, header=header)

    def add(self, hdu_or_table, sub_pixel=False, order="bilinear"):
        """
        Add a projection of an image or table files to the canvas

        .. note::
          If a Table is provided, it must include the following columns:
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
            Default is False. Dictates if point files should be projected with
            sub-pixel shifts or not. Accounting for sub-pixel shifts is approx.
            5x slower.

        """

        if isinstance(hdu_or_table, Table):
            self.hdu.header["COMMENT"] = "Adding files from table"
            self.hdu = add_table_to_imagehdu(hdu_or_table, self.hdu,
                                             sub_pixel=sub_pixel)
        elif isinstance(hdu_or_table, fits.ImageHDU):
            self.hdu.header["COMMENT"] = "Adding files from table"
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


