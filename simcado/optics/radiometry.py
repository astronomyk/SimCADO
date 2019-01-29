import os
import warnings
from collections import OrderedDict
from copy import deepcopy

import numpy as np
from astropy.table import Table, Row, vstack
from astropy.io import ascii as ioascii
from astropy import units as u
from synphot import SpectralElement, SourceSpectrum
from synphot import Empirical1D

from .surface import SpectralSurface, quantify
from ..server.database import change_table_entry
from .. import utils


class RadiometryTable:

    def __init__(self, tables=(), **kwargs):
        self.meta = {}
        self.meta.update(kwargs)

        self.table = None
        self.surfaces = OrderedDict({})

        if len(tables) > 0:
            self.add_surface_list(tables)

    def add_surface_list(self, surface_lists, prepend=False):
        self.table = combine_tables(surface_lists, self.table, prepend=prepend)

        r_name = real_colname("name", self.table.colnames)
        for row in self.table:
            if row[r_name] not in self.surfaces:
                surf = make_surface_from_row(row)
                self.add_surface(surf, row[r_name], position=-1,
                                 add_to_table=False)

    def add_surface(self, surface, name, position=-1, add_to_table=True):
        if self.table is None:
            raise ValueError("Cannot add surface without <self>.table template."
                             "Please add an empty table to define column names")

        self.surfaces = add_surface_to_dict(self.surfaces, surface,
                                            name  , position)
        if add_to_table:
            self.table = add_surface_to_table(self.table, surface,
                                              name, position)

    def get_throughput(self, start=0, end=None, rows=None):
        if self.table is None:
            return None

        if end is None:
            end = len(self.table)
        elif end < 0:
            end += len(self.table)
        if rows is None:
            rows = np.arange(start, end)

        return combine_throughputs(self.table, self.surfaces, rows)

    def get_emission(self, etendue, start=0, end=None, rows=None):
        if self.table is None:
            return None

        if end is None:
            end = len(self.table)
        elif end < 0:
            end += len(self.table)
        if rows is None:
            rows = np.arange(start, end)

        return combine_emissions(self.table, self.surfaces, rows, etendue)

    @property
    def emission(self):
        if "etendue" not in self.meta["etendue"]:
            raise ValueError("self.meta['etendue'] must be set")
        etendue = quantify(self.meta["etendue"], u.Unit("m2 arcsec2"))

        return self.get_emission(etendue)

    @property
    def throughput(self):
        return self.get_throughput()

    def plot(self, what="all", rows=None):
        raise NotImplemented

    def __getitem__(self, item):
        return self.surfaces[item]

    def __repr__(self):
        print(self.table)


def combine_emissions(tbl, surfaces, row_indexes, etendue):
    if len(tbl) == 0:
        return None

    etendue = quantify(etendue, "m2 arcsec2")

    r_name = real_colname("name", tbl.colnames)
    r_action = real_colname("action", tbl.colnames)

    emission = None
    for ii, row_num in enumerate(row_indexes):
        row = tbl[row_num]
        surf = surfaces[row[r_name]]
        action_attr = row[r_action]

        if isinstance(surf, SpectralSurface):
            surf_throughput = getattr(surf, action_attr)

            surf_emission = surf.emission
            surf_eff_area = surf.area * np.cos(surf.mirror_angle)
            surf_eff_solid_angle = (etendue / surf_eff_area).to(u.arcsec**2)
            surf_emission *= surf_eff_solid_angle.value

            surf_emission.meta["solid_angle"] = None
            surf_emission.meta["history"] += ["Etendue scale factor applied. "
                                              "Effective pixel solid angle for "
                                              "surface is {}"
                                              "".format(surf_eff_solid_angle)]

            if ii == 0:
                emission = deepcopy(surf_emission)
            else:
                emission = emission * surf_throughput
                emission = emission + surf_emission

    return emission


def combine_throughputs(tbl, surfaces, rows_indexes):
    if len(tbl) == 0:
        return None

    r_name = real_colname("name", tbl.colnames)
    r_action = real_colname("action", tbl.colnames)

    throughput = None
    for ii, row_num in enumerate(rows_indexes):

        row = tbl[row_num]
        surf = surfaces[row[r_name]]
        action_attr = row[r_action]

        if isinstance(surf, SpectralSurface):
            surf_throughput = getattr(surf, action_attr)

            if ii == 0:
                throughput = deepcopy(surf_throughput)
            else:
                throughput = throughput * surf_throughput

    return throughput


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
            if isinstance(surf_val, u.Quantity):
                surf_val = surf_val.value
            tbl = change_table_entry(tbl, colname, surf_val, position=position)

    colname = real_colname("name", tbl.colnames)
    tbl = change_table_entry(tbl, colname, name, position=position)

    return tbl


def add_surface_to_dict(dic, surf, name, position=0):
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
