import numpy as np
from astropy import units as u

from ... import utils
from .effects import Effect
from ..image_plane_utils import _header_from_list_of_xy


class DetectorList(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)

    @property
    def image_plane_header(self):
        tbl = self.get_data()
        pixel_size = np.min(utils.quantity_from_table("pixsize", tbl, u.mm))
        x_unit = utils.unit_from_table("x_cen", tbl, u.mm)
        y_unit = utils.unit_from_table("y_cen", tbl, u.mm)

        x_det_min = np.min(tbl["x_cen"] - tbl["xhw"]) * x_unit
        x_det_max = np.max(tbl["x_cen"] + tbl["xhw"]) * x_unit
        y_det_min = np.min(tbl["y_cen"] - tbl["yhw"]) * y_unit
        y_det_max = np.max(tbl["y_cen"] + tbl["yhw"]) * y_unit

        x_det = [x_det_min.to(u.mm).value, x_det_max.to(u.mm).value]
        y_det = [y_det_min.to(u.mm).value, y_det_max.to(u.mm).value]

        pixel_size = pixel_size.to(u.mm).value
        hdr = _header_from_list_of_xy(x_det, y_det, pixel_size, "D")

        return hdr

    def detector_headers(self, ids=None):
        """Return the header to describe each individual detector"""
        raise NotImplementedError
