"""
Put functions which connect the IPDB to SimCADO here
Possibly also add the skycalc_cli interface here

Write functions to:
0. Connect to the server
1a. Read which packages are available
1b. Look at which packages are available locally
2. Display a list of which packages are available
3. Download a package
4. Unpack into it's own folder

"""

import os
import requests
import datetime as dt
from copy import deepcopy

from numpy import where as npwhere
from numpy import array as nparray
from astropy.table import Column, Row, vstack
from astropy.io import ascii as ioascii

from .. import default_keywords as dkeys


LOCAL_DB_PATH   = os.path.join(dkeys.PKG_DIR,
                               dkeys.INST_PKG_LOCAL_PATH,
                               dkeys.INST_PKG_LOCAL_DB_NAME)
SERVER_DB_PATH  = os.path.join(dkeys.INST_PKG_SERVER_PATH,
                               dkeys.INST_PKG_SERVER_DB_NAME)



def get_local_packages(path=None):
    if path is None:
        path = LOCAL_DB_PATH

    if not os.path.exists(path):
        raise ValueError(path + " doesn't exist")

    # If this throws an error, it's because the DB file is not formated correcty
    # The column row names should NOT be commented out
    local_table = ioascii.read(path, format="basic")
    local_table = rename_table_colnames(local_table)

    col = Column(name="origin", data=["local"] * len(local_table))
    local_table.add_column(col)

    return local_table


def get_server_text(path=None):
    if path is None:
        path = SERVER_DB_PATH

    server_db_text = requests.get(path).text
    server_db_text = server_db_text.replace("\r", "")

    return server_db_text


def rename_table_colnames(svr_table):
    if svr_table.colnames[0] != "col0":
        return svr_table

    col_names = svr_table.meta["comments"][-1].replace("#", "").split()

    if len(col_names) != len(svr_table.colnames):
        raise ValueError("Somethings up with the server database header names")

    for i, col_name in enumerate(col_names):
        svr_table[svr_table.colnames[i]].name = col_name

    return svr_table


def get_server_packages(path=None):
    svr_text = get_server_text(path)
    svr_table = ioascii.read(svr_text)

    svr_table = rename_table_colnames(svr_table)

    col = Column(name="origin", data=["server"] * len(svr_table))
    svr_table.add_column(col)

    return svr_table


def list_packages(local_path=None, server_url=None, return_table=False):
    """
    Prints to screen (returns) a list of all available packages
        
    Parameters
    ----------
    local_path : str, optional
    server_url : str, optional
    return_table : bool, optional

    Returns
    -------
    if return_table == True
        astropy.Table

    """

    local_table = get_local_packages(local_path)
    svr_table = get_server_packages(server_url)
    all_table = vstack(local_table, svr_table)

    print("\nPackages saved offline\n======================")
    print(local_table)
    print("\nPackages on the server\n======================")
    print(svr_table)

    if return_table:
        return all_table


def check_package_exists(pkg_name, svr_path=None):
    if svr_path is None:
        svr_path = SERVER_DB_PATH

    svr_base_url = os.path.dirname(svr_path)
    svr_table = get_server_packages(svr_path)
    pkg_entry = get_package_table_entry(pkg_name, svr_table)
    url = os.path.join(svr_base_url, pkg_entry["path"]).replace("\\", "/")

    if not requests.get(url).ok:
        raise ValueError(url + " doesn't return ALL_GOOD (200)")

    return True


def get_server_package_path(pkg_name, svr_table):
    if pkg_name not in svr_table["name"]:
        raise ValueError(pkg_name + " not in server table")

    svr_dict = {n: p for n, p in zip(svr_table["name"], svr_table["path"])}
    pkg_path = svr_dict[pkg_name]

    return pkg_path


def get_package_table_entry(pkg_name, svr_table, multiple_entries="newest"):

    pkg_index = npwhere(svr_table["name"] == pkg_name)[0]
    if pkg_index.shape[0] == 0:
        raise ValueError(pkg_name+" was not found in the table")

    # ::todo implement multiple entry handling
    pkg_entry = svr_table[pkg_index[0]]

    return pkg_entry


def download_package(pkg_name, save_dir=None, svr_db_path=None):
    """
    Download a package from the server

    Parameters
    ----------
    pkg_name : str
        Name of the package to download

    save_dir : str, optional
        Where to save the package on the local disk. Default INST_PKG_LOCAL_PATH

    svr_db_path : str, optional
        URL to the server

    """

    if save_dir is None:
        save_dir = os.path.join(dkeys.PKG_DIR,
                                dkeys.INST_PKG_LOCAL_PATH)

    if svr_db_path is None:
        svr_path = SERVER_DB_PATH

    if not check_package_exists(pkg_name):
        raise ValueError("Package doesn't exist: " + pkg_name)

    svr_base_url = os.path.dirname(svr_db_path)
    svr_table = get_server_packages(svr_db_path)
    pkg_entry = get_package_table_entry(pkg_name, svr_table)
    pkg_url = os.path.join(svr_base_url, pkg_entry["path"]).replace("\\", "/")

    from ..utils.utils import download_file
    local_filename = download_file(pkg_url, save_dir)
    print("Saved {} in {}".format(pkg_name, local_filename))

    pkg_entry_local = deepcopy(pkg_entry)
    pkg_entry_local["path"] = local_filename


def add_pkg_to_local_db(new_pkg_entry, local_table=None):

    if local_table is None:
        local_table = get_local_packages(LOCAL_DB_PATH)
    if type(new_pkg_entry) != Row:
        raise ValueError("pkg_entry must be a astropy.table.Row object")

    print(new_pkg_entry, local_table)

    if new_pkg_entry["name"] in local_table["name"]:
        ii = npwhere(local_table["name"] == new_pkg_entry["name"])[0][0]
        fmt = '%Y-%m-%d'
        local_date = dt.datetime.strptime(local_table[ii]["date_modified"], fmt)
        new_date   = dt.datetime.strptime(new_pkg_entry["date_modified"], fmt)

        if new_date >= local_date:
            local_table[ii]["name"] = local_table[ii]["name"] + "_" + \
                                      local_table[ii]["date_modified"]
            print(local_table[ii]["name"])
        else:


            new_pkg_entry["name"] = new_pkg_entry["name"] + "_" + \
                                    new_pkg_entry["date_modified"]
            print(new_pkg_entry[ii]["name"])
        #print(ii, local_table[ii]["name"], new_pkg_entry["name"],)

    local_table.add_row(new_pkg_entry)
    print(local_table)

    return local_table


def change_table_entry(tbl, col_name, old_val, new_val):

    offending_col = list(tbl[col_name].data)

    for ii in npwhere(old_val in offending_col)[0]:
        offending_col[ii] = new_val

    fixed_col = Column(name=col_name, data=offending_col)

    ii = npwhere(nparray(tbl.colnames) == col_name)[0][0]
    tbl.remove_column(col_name)
    tbl.add_column(fixed_col, index=ii)

    return tbl, fixed_col, offending_col, new_val, ii

















# check_local_directory_exists():






# def add_pkg_to_table(pkg_name, local_filename, svr_table):
#
#
