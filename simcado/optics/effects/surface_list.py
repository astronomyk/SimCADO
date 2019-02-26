import numpy as np

from ... import rc

from .effects import Effect, TERCurve
from ..radiometry import RadiometryTable


class SurfaceList(Effect):
    def __init__(self, **kwargs):
        super(SurfaceList, self).__init__(**kwargs)

        self.meta["SIM_MIN_THRESHOLD"] = rc.__rc__["SIM_MIN_THRESHOLD"]

        self.radiometry_table = RadiometryTable()
        self.radiometry_table.meta.update(self.meta)

        data = self.get_data()
        if data is not None:
            self.radiometry_table.add_surface_list(data)

    def apply_to(self):
        pass

    def fov_grid(self, header, waverange):
        wave = np.linspace(min(waverange), max(waverange), 100)
        throughput = self.get_throughput()(wave)
        valid_waves = np.where(throughput > self.meta["SIM_MIN_THRESHOLD"])[0]
        if len(valid_waves) > 0:
            wave_edges = [min(wave[valid_waves]), max(wave[valid_waves])]
        else:
            raise ValueError("No transmission found above the threshold {} in "
                             "this wavelength range: {}"
                             "".format(self.meta["SIM_MIN_THRESHOLD"],
                                       waverange))

        return {"coords": None, "wavelengths": wave_edges}

    def add_surface(self, surface, name, position=-1, add_to_table=True):
        if isinstance(surface, TERCurve):
            surface = surface.surface
        self.radiometry_table.add_surface(surface, name, position, add_to_table)

    def add_surface_list(self, surface_list, prepend=False):
        if isinstance(surface_list, SurfaceList):
            surface_list = surface_list.radiometry_table.table
        elif isinstance(surface_list, RadiometryTable):
            surface_list = surface_list.table

        self.radiometry_table.add_surface_list(surface_list, prepend)

    def get_emission(self, **kwargs):
        if "etendue" in kwargs:
            etendue = kwargs["etendue"]
        elif "etendue" in self.meta:
            etendue = self.meta["etendue"]
        elif "etendue" in self.radiometry_table.meta:
            etendue = self.radiometry_table.meta["etendue"]
        else:
            raise ValueError("etendue must be given in kwargs or .meta")

        return self.radiometry_table.get_emission(etendue=etendue, **kwargs)

    def get_throughput(self, **kwargs):
        return self.radiometry_table.get_throughput(**kwargs)
