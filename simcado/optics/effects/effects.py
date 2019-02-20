from simcado.optics.fov import FieldOfView
from ..data_container import DataContainer


class Effect(DataContainer):

    def __init__(self, **kwargs):
        super(Effect, self).__init__(**kwargs)

    def apply_to(self, fov):
        if not isinstance(fov, FieldOfView):
            raise ValueError("fov must be a FieldOfView object: {}"
                             "".format(type(fov)))

    def fov_grid(self, header, waverange):
        return {"coords": None, "wavelengths": None}

    def update(self, **kwargs):
        pass


class SurfaceList(Effect):
    def __init__(self, **kwargs):
        super(SurfaceList, self).__init__(**kwargs)


class ApertureList(Effect):
    def __init__(self, **kwargs):
        super(ApertureList, self).__init__(**kwargs)
