from simcado.optics.fov import FieldOfView
from ..data_container import DataContainer
from ..surface import SpectralSurface


class Effect(DataContainer):

    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)

    def apply_to(self, fov):
        if not isinstance(fov, FieldOfView):
            raise ValueError("fov must be a FieldOfView object: {}"
                             "".format(type(fov)))

    def fov_grid(self, header, waverange, **kwargs):
        return {"coords": None, "wavelengths": None}

    def update(self, **kwargs):
        pass

    def __repr__(self):
        if "name" not in self.meta:
            self.meta["name"] = "<unknown name>"
        return '{}: "{}"'.format(type(self).__name__, self.meta["name"])


class TERCurve(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)

        self.surface = SpectralSurface()
        self.surface.meta.update(self.meta)
        data = self.get_data()
        if data is not None:
            self.surface.table = data


class Shift3D(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)


class AtmosphericDispersion(Shift3D):
    def __init__(self, **kwargs):
        super(Shift3D, self).__init__(**kwargs)

    def fov_grid(self, header, waverange, **kwargs):
        pass


class AtmosphericDispersionCorrector(Shift3D):
    def __init__(self, **kwargs):
        super(Shift3D, self).__init__(**kwargs)


class ApertureList(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)


class ApertureMask(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)
        self.header = None


class TraceList(Effect):
    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)




