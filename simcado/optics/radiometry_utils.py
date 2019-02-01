from collections import OrderedDict
from copy import deepcopy

import numpy as np
from astropy import units as u
from astropy.io import ascii as ioascii
from astropy.table import Table, vstack

from .surface import SpectralSurface
from ..utils import real_colname, insert_into_ordereddict, quantify, \
    change_table_entry, convert_table_comments_to_dict


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
    meta_dict = convert_table_comments_to_dict(tbl)
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