import os
import warnings
from collections import OrderedDict

from astropy.table import Table, Row, vstack
from astropy.io import ascii as ioascii

from .surface import SpectralSurface
from ..server.database import change_table_entry
from .. import utils


class RadiometryTable:

    def __init__(self, tables=(), surfaces=(), **kwargs):
        self.meta = {}
        self.meta.update(kwargs)

        self.table = None
        self.surfaces = OrderedDict({})

        if len(tables) > 0:
            self.add_surface_list(tables)

    def add_surface(self, surface, name, position):
        self.surfaces = add_surface_to_dict(self.surfaces, surface,
                                            name, position)
        add_surface_to_table(self.surfaces, surface,
                                            name, position)

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


def add_surface_to_table(tbl, surf, name, position):
    tbl.insert_row(position)
    for colname in tbl.colnames:
        surf_col = real_colname(colname, surf.meta)
        if surf_col:
            surf_val = surf.meta[surf_col]
            tbl = change_table_entry(tbl, colname, surf_val, position=position)

    tbl_name_col = real_colname("name", tbl.colnames)
    tbl[tbl_name_col][position] = name

    return tbl


def add_surface_to_dict(dic, surf, name, position=-1):
    new_entry = OrderedDict({name : surf})
    dic = insert_into_ordereddict(dic, new_entry, position)

    return dic


def make_surface_dict_from_table(tbl):
    names = tbl[real_colname("name", tbl.colnames)]
    surf_dict = OrderedDict({})
    for ii in range(len(tbl)):
        surf_dict[names[ii]] = make_surface_from_row(tbl[ii], **tbl.meta)

    return surf_dict


def make_surface_from_row(row, **kwargs):
    row_dict = {colname.lower() : row[colname] for colname in row.colnames}
    kwargs.update(row_dict)
    surface = SpectralSurface(**kwargs)

    return surface


def real_colname(name, colnames):
    names = [name.lower(), name.upper(), name[0].upper() + name[1:].lower()]
    real_name = [name for name in names if name in colnames]
    if len(real_name) == 0:
        real_name = None
    else:
        real_name = real_name[0]

    return real_name


def insert_into_ordereddict(dic, new_entry, pos):
    if isinstance(new_entry, dict):
        new_entry = [[key, val] for key, val in new_entry.items()]
    elif isinstance(new_entry, (list, tuple)) and \
            not isinstance(new_entry[0], (list, tuple)):
        new_entry = [new_entry]

    if pos < 0:
        pos += len(dic) + len(new_entry)

    new_dic = list(OrderedDict(dic).items())
    new_dic = new_dic[:pos] + new_entry + new_dic[pos:]
    new_dic = OrderedDict(new_dic)

    return new_dic


def empty_type(x):
    type_dict = {int: 0, float: 0., bool: False, str: " ",
                 list: [], tuple: (), dict: {}}
    if "<U" in str(x):
        x = str

    return type_dict[x]

