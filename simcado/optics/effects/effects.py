from astropy.io import fits

from ...source.source2 import Source
from ..fov import FieldOfView
from ..data_container import DataContainer
from ..surface import SpectralSurface
from ..image_plane_utils import calc_footprint


class Effect(DataContainer):

    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)

    def apply_to(self, obj):
        if obj is None:
            raise ValueError("fov was None. Oops.")

        if not isinstance(obj, (Source, FieldOfView)):
            raise ValueError("fov must be a FieldOfView (or Source) object: {}"
                             "".format(type(obj)))

        return obj

    def fov_grid(self, header=None, waverange=None, **kwargs):
        return {"edges": None, "wavelengths": None}

    def update(self, **kwargs):
        pass

    def __repr__(self):
        if "name" not in self.meta:
            self.meta["name"] = "<unknown name>"
        return '{}: "{}"'.format(type(self).__name__, self.meta["name"])


class TERCurve(Effect):
    def __init__(self, **kwargs):
        super(TERCurve, self).__init__(**kwargs)

        self.surface = SpectralSurface()
        self.surface.meta.update(self.meta)
        data = self.get_data()
        if data is not None:
            self.surface.table = data


class AtmosphericTERCurve(TERCurve):
    def __init__(self, **kwargs):
        super(AtmosphericTERCurve, self).__init__(**kwargs)


class Shift3D(Effect):
    def __init__(self, **kwargs):
        super(Shift3D, self).__init__(**kwargs)

    def fov_grid(self, header=None, waverange=None, **kwargs):
        dic = {"wavelengths": waverange, "x_shifts": [0, 0], "y_shifts": [0, 0]}
        return dic


class AtmosphericDispersion(Shift3D):
    def __init__(self, **kwargs):
        super(AtmosphericDispersion, self).__init__(**kwargs)


class AtmosphericDispersionCorrector(Shift3D):
    def __init__(self, **kwargs):
        super(AtmosphericDispersionCorrector, self).__init__(**kwargs)


class ApertureList(Effect):
    def __init__(self, **kwargs):
        super(ApertureList, self).__init__(**kwargs)


class ApertureMask(Effect):
    def __init__(self, **kwargs):
        super(ApertureMask, self).__init__(**kwargs)

    def fov_grid(self, header=None, waverange=None, **kwargs):
        x_sky, y_sky = calc_footprint(self.header)
        return {"wavelengths": None, "edges": [x_sky, y_sky]}

    @property
    def header(self):
        return fits.Header()


class TraceList(Effect):
    def __init__(self, **kwargs):
        super(TraceList, self).__init__(**kwargs)




