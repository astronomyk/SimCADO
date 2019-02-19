from ..data_container import DataContainer


class Effect(DataContainer):

    def __init__(self, **kwargs):
        self.meta = {}
        self.meta.update(kwargs)
        super(Effect, self).__init__(**kwargs)

    def apply_to(self, fov):
        return fov

    def fov_grid(self, header, waverange):
        coords = None
        wavelengths = None

        return coords, wavelengths
