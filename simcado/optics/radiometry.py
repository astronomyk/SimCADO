from astropy import table

from .surface import SpectralSurface


class RadiometryTable:

    def __init__(self, tables=()):
        self.tbl = table.Table()
        self.surfaces = {}

    def add_surface(self, ):
        pass

    def add_surface_list(self):
        pass

    def get_transmission(self, rows):
        pass

    def get_emission(self, rows=None, ):
        pass

    def list(self, what="all"):
        pass

    def plot(self, what="all", rows=None):
        pass

    def __getitem__(self, item):
        pass

