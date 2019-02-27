from astropy import units as u
from astropy.table import Table

from simcado.optics.effects import SurfaceList, TERCurve


def _surf_list():
    kwargs = {"etendue": 5776 * u.m ** 2 * u.mas ** 2,
              "filename": "EC_mirrors_MICADO_Wide.tbl",
              "name": "MICADO Mirror List"}
    return SurfaceList(**kwargs)


def _surf_list_empty():
    tbl = Table(names=["Name", "Outer", "Inner", "Angle",
                       "Temp", "Action", "Filename"],
                meta={"outer_unit": "m", "inner_unit": "m", "angle_unit": "deg",
                      "temp_unit": "deg_C"})
    kwargs = {"etendue": 5776 * u.m ** 2 * u.mas ** 2,
              "table": tbl,
              "name": "Empty Surface List"}
    return SurfaceList(**kwargs)


def _filter_surface():
    kwargs = {"filename": "TC_filter_Ks.dat",
              "name": "filter",
              "action": "transmission",
              "outer": 0.1,
              "temp": 0}
    return TERCurve(**kwargs)
