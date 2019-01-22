import os
import warnings
from _collections import OrderedDict

from astropy.table import Table, Row, vstack
from astropy.io import ascii as ioascii

from .surface import SpectralSurface
from .. import utils


class RadiometryTable:

    def __init__(self, tables=(), **kwargs):
        self.meta = {}
        self.meta.update(kwargs)

        self.table = None
        self.surfaces = {}

        if len(tables) > 0:
            self.add_surface_list(tables)

    def add_surface(self, surface):
        pass

    def add_surface_list(self, surface_lists, prepend=False):
        self.table = combine_tables(surface_lists, self.table, prepend=prepend)

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


def combine_tables(new_tables, old_table=None, prepend=False):
    if isinstance(old_table, str):
        old_table = string_to_table(old_table)

    if isinstance(new_tables, (str, Table)):
        new_tables = [new_tables]

    for new_table in new_tables:
        if isinstance(new_table, str):
            new_table = string_to_table(new_table)
        if old_table is None:
            old_table = new_table
        else:
            if prepend:
                old_table = vstack([new_table, old_table])
            else:
                old_table = vstack([old_table, new_table])

    return old_table


def string_to_table(tbl):
    tbl = ioascii.read(tbl, fast_reader=False)
    meta_dict = utils.convert_table_comments_to_dict(tbl)
    tbl.meta.update(meta_dict)

    return tbl


def add_surface_to_table():
    pass


def make_surface_dict_from_table(tbl):
    names = tbl[real_colname("name", tbl.colnames)]
    surf_dict = OrderedDict({})
    for ii in range(len(tbl)):
        surf_dict[names[ii]] = make_surface_from_row(tbl[ii], **tbl.meta)

    return surf_dict


def make_surface_from_row(row, **kwargs):
    row_dict = {colname.lower() : row[colname] for colname in row.colnames}
    row_dict.update(kwargs)
    surface = SpectralSurface(**row_dict)

    return surface


def real_colname(name, colnames):
    names = [name.lower(), name.upper(), name[0].upper() + name[1:].lower()]
    real_name = [name for name in names if name in colnames][0]

    return real_name

