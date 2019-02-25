import numpy as np
from astropy import units as u

from ... import utils
from .effects import Effect
from ..image_plane_utils import _header_from_list_of_xy


class DetectorList(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)

    def image_plane_header(self, pixel_scale):
        """
        Make the header to describe the extent of the image plane

        Parameters
        ----------
        pixel_scale : u.Quantity
            [arcsec]

        Returns
        -------
        hdr : fits.Header

        """

        tbl = self.get_data()

        pixel_scale = utils.quantify(pixel_scale, u.arcsec)       # arcsec pix-1
        pixel_size = utils.quantity_from_table("pixsize", tbl, u.mm)  # mm pix-1
        plate_scale = pixel_scale / pixel_size                     # arcsec mm-1

        x_unit = utils.unit_from_table("x_cen", tbl, u.mm)
        y_unit = utils.unit_from_table("y_cen", tbl, u.mm)

        x_sky_min = np.min((tbl["x_cen"] - tbl["xhw"]) * plate_scale) * x_unit
        x_sky_max = np.max((tbl["x_cen"] + tbl["xhw"]) * plate_scale) * x_unit
        y_sky_min = np.min((tbl["y_cen"] - tbl["yhw"]) * plate_scale) * y_unit
        y_sky_max = np.max((tbl["y_cen"] + tbl["yhw"]) * plate_scale) * y_unit

        x_sky = [x_sky_min.to(u.deg).value, x_sky_max.to(u.deg).value]
        y_sky = [y_sky_min.to(u.deg).value, y_sky_max.to(u.deg).value]

        pixel_scale = pixel_scale.to(u.deg).value
        hdr = _header_from_list_of_xy(x_sky, y_sky, pixel_scale)

        return hdr

    def detector_headers(self, ids=None):
        """Return the header to describe each individual detector"""
        raise NotImplementedError
