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

from numpy import where as npwhere
from numpy import array as nparray
from astropy.table import Table, Column, Row, vstack
from astropy.io import ascii as ioascii

from .. import __rc__
rc = __rc__


def set_local_path_names(path):
    local_inst_db = os.path.join(path, rc["FILE_INST_PKG_LOCAL_DB_NAME"])
    local_psf_db = os.path.join(path, rc["FILE_PSF_LOCAL_DB_NAME"])
    local_src_db = os.path.join(path, rc["FILE_SRC_PKG_LOCAL_DB_NAME"])

    return local_inst_db, local_psf_db, local_src_db


_local_paths = set_local_path_names(rc["FILE_LOCAL_DOWNLOADS_PATH"])
LOCAL_INST_DB = _local_paths[0]
LOCAL_PSF_DB  = _local_paths[1]
LOCAL_SRC_DB  = _local_paths[2]

SVR_INST_DB = rc["FILE_SERVER_BASE_URL"] + rc["FILE_INST_PKG_SERVER_DB_NAME"]
SVR_PSF_DB  = rc["FILE_SERVER_BASE_URL"] + rc["FILE_PSF_SERVER_DB_NAME"]
SVR_SRC_DB  = rc["FILE_SERVER_BASE_URL"] + rc["FILE_SRC_PKG_SERVER_DB_NAME"]

_local_db_dict = {"inst" : LOCAL_INST_DB,
                  "psf"  : LOCAL_PSF_DB,
                  "src"  : LOCAL_SRC_DB}

_svr_db_dict = {"inst" : SVR_INST_DB,
                "psf"  : SVR_PSF_DB,
                "src"  : SVR_SRC_DB}


LOCAL_DB_HEADER_PATTERN = """# Date-created : {}
# Date-modified : {}
# Description : Packages containing {} specific data files
#
name   author   date_added  date_modified   path
.      .        .           .               .     
"""


def set_up_local_package_directory(dirname=None, overwrite=False):

    if dirname is None:
        dirname = rc["FILE_LOCAL_DOWNLOADS_PATH"]

    # set up downloads directory and sub directories
    for dname in [dirname,
                  os.path.join(dirname, rc["FILE_SCOPE_PKG_LOCAL_PATH"]),
                  os.path.join(dirname, rc["FILE_INST_PKG_LOCAL_PATH"]),
                  os.path.join(dirname, rc["FILE_PSF_LOCAL_PATH"]),
                  os.path.join(dirname, rc["FILE_SRC_PKG_LOCAL_PATH"])]:
        if not os.path.exists(dname):
            os.makedirs(dname)
        elif not overwrite:
            print("{} already exists. If you would like to overwrite it, set"
                  "overwrite=True".format(dname))

    # make db files
    now = dt.datetime.now().strftime('%Y-%m-%d')
    local_paths = set_local_path_names(dirname)
    package_type = ["instrument/telescope", "PSF", "Source object"]

    for loc_db, pkg_type in zip(local_paths, package_type):
        if not os.path.exists(loc_db) or overwrite:
            with open(loc_db, "w") as new_db:
                hdr_text = LOCAL_DB_HEADER_PATTERN.format(now, now, pkg_type)
                new_db.write(hdr_text)
        elif not overwrite:
            print("{} already exists. If you would like to overwrite it, set"
                  "'overwrite=True'".format(loc_db))


def get_local_packages(path=None):
    """
    Returns a table of packages from a database file

    Parameters
    ----------
    path : str
        Path to package database file. Use the following local variables:
        `simcado.server.LOCAL_INST_DB`
        `simcado.server.LOCAL_PSF_DB`
        `simcado.server.LOCAL_SRC_DB`

    Returns
    -------
    local_table : `astropy.Table`

    """

    if path is None:
        path = LOCAL_INST_DB

    if not os.path.exists(path):
        raise ValueError(path + " doesn't exist")

    # If this throws an error, it's because the DB file is not formatted
    # correctly
    # The column row names should NOT be commented out
    local_table = ioascii.read(path, format="basic")
    local_table = rename_table_colnames(local_table)
    local_table = remove_dot_row(local_table)

    # make sure all columns are strings
    # for col in local_table.colnames:
    #    local_table[col] = local_table[col].astype("str")

    return local_table


def remove_dot_row(tbl):
    if "." in tbl["name"]:
        ii = npwhere(tbl["name"] == ".")[0][0]
        tbl.remove_row(ii)
    return tbl


def get_server_text(path=None):
    if path is None:
        path = SVR_INST_DB

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

    return svr_table


def list_packages(local_path=None, server_url=None, return_table=False,
                  msg="Packages"):
    """
    Prints to screen (returns) a list of all available packages
        
    Parameters
    ----------
    local_path : str, optional
        Default `simcado.server.LOCAL_PSF_DB`

    server_url : str, optional
        Default `simcado.server.SVR_PSF_DB`

    return_table : bool, optional

    msg : str


    Returns
    -------
    all_table : `astropy.Table`
        Only if `return_table` is `True`

    """

    local_table = get_local_packages(local_path)
    svr_table = get_server_packages(server_url)
    all_table = vstack(local_table, svr_table)

    print("\n{} saved offline\n".format(msg) + "="*(len(msg)+14))
    print(local_table)
    print("\n{} on the server\n".format(msg) + "="*(len(msg)+14))
    print(svr_table)

    if return_table:
        return all_table


def list_instruments(local_path=None, server_url=None, return_table=False):
    """
    Prints to screen (returns) a list of all available instrument packages

    By default `list_instruments` looks in

    * `simcado.server.LOCAL_INST_DB`
    * `simcado.server.SVR_INST_DB`

    See Also
    --------
    :func:`.list_packages`

    """

    if local_path is None:
        local_path = LOCAL_INST_DB
    if server_url is None:
        server_url = SVR_INST_DB

    return list_packages(local_path, server_url, return_table,
                         "Instrument and Telescope packages")


def list_psfs(local_path=None, server_url=None, return_table=False):
    """
    Prints to screen (returns) a list of all available PSF files

    By default `list_psfs` looks in

    * `simcado.server.LOCAL_PSF_DB`
    * `simcado.server.SVR_PSF_DB`

    See Also
    --------
    :func:`.list_packages`

    """

    if local_path is None:
        local_path = LOCAL_PSF_DB
    if server_url is None:
        server_url = SVR_PSF_DB

    return list_packages(local_path, server_url, return_table,
                         "PSF files")


def list_source_pkgs(local_path=None, server_url=None, return_table=False):
    """
    Prints to screen (returns) a list of all available source packages

    By default `list_source_pkgs` looks in

    * `simcado.server.LOCAL_SRC_DB`
    * `simcado.server.SVR_SRC_DB`

    See Also
    --------
    :func:`.list_packages`

    """

    if local_path is None:
        local_path = LOCAL_SRC_DB
    if server_url is None:
        server_url = SVR_SRC_DB

    return list_packages(local_path, server_url, return_table,
                         "Source packages")


def list_all():
    list_instruments()
    list_psfs()
    list_source_pkgs()


def check_package_exists(pkg_name, svr_path=None):
    if svr_path is None:
        svr_path = SVR_INST_DB

    svr_base_url = os.path.dirname(svr_path)
    svr_table = get_server_packages(svr_path)
    pkg_entry = get_package_table_entry(pkg_name, svr_table)

    if pkg_entry is not None:
        url = os.path.join(svr_base_url, pkg_entry["path"]).replace("\\", "/")

        if not requests.get(url).ok:
            raise ValueError(url + " doesn't return ALL_GOOD (200)")

        return_val = True
    else:
        return_val = False

    return return_val


def get_server_package_path(pkg_name, svr_table):

    pkg_path = None
    if pkg_name in svr_table["name"]:
        svr_dict = {n: p for n, p in zip(svr_table["name"], svr_table["path"])}
        pkg_path = svr_dict[pkg_name]

    return pkg_path


def get_package_table_entry(pkg_name, db_table):
    # ::todo implement multiple entry handling
    if pkg_name in db_table["name"]:
        pkg_index = npwhere(db_table["name"] == pkg_name)[0]
        pkg_entry = db_table[pkg_index[0]]
    else:
        pkg_entry = None

    return pkg_entry


def determine_type_of_package(svr_db_filename):

    pkg_type = None
    if "inst" in svr_db_filename.lower():
        pkg_type = "inst"
    elif "psf" in svr_db_filename.lower():
        pkg_type = "psf"
    elif "source" in svr_db_filename.lower():
        pkg_type = "src"

    return pkg_type


def download_package(pkg_name, save_dir=None, server_dbs=None):
    """
    Download a package from the server

    Parameters
    ----------
    pkg_name : str
        Name of the package to download

    save_dir : str, optional
        Where to save the package on the local disk. Default INST_PKG_LOCAL_PATH

    server_dbs : str, optional
        URL to the server

    """

    pkg_entry, svr_db = find_package_on_server(pkg_name,
                                               server_dbs=server_dbs,
                                               return_db_filename=True)

    if pkg_entry is None:
        raise ValueError("{} wasn't found on the server".format(pkg_name))

    pkg_url  = rc["FILE_SERVER_BASE_URL"] + pkg_entry["path"]
    pkg_type = determine_type_of_package(svr_db)

    if not check_package_exists(pkg_name, _svr_db_dict[pkg_type]):
        raise ValueError("Package is missing: " + pkg_name)

    if save_dir is None:
        stem = os.path.dirname(pkg_entry["path"])
        save_dir = os.path.join(rc["FILE_LOCAL_DOWNLOADS_PATH"], stem)

    from ..utils import download_file
    local_filename = download_file(pkg_url, save_dir)
    print("Saved {} in {}".format(pkg_name, local_filename))

    local_db_path = _local_db_dict[pkg_type]
    new_local_tbl = add_pkg_to_local_db(pkg_entry, local_db_path)

    write_table_to_disk(new_local_tbl, local_db_path)


def write_table_to_disk(tbl, path):

    tbl.write(path, format="ascii.fixed_width", overwrite=True, delimiter="")


def find_package_on_server(pkg_name, server_dbs=None, return_db_filename=False):

    if server_dbs is None:
        server_dbs = [SVR_INST_DB, SVR_PSF_DB, SVR_SRC_DB]

    pkg_entry, svr_db = None, None
    for svr_db in server_dbs:
        svr_table = get_server_packages(svr_db)
        pkg_entry = get_package_table_entry(pkg_name, svr_table)
        if pkg_entry is not None:
            break

    if return_db_filename:
        return pkg_entry, os.path.basename(svr_db)
    else:
        return pkg_entry


def add_pkg_to_local_db(new_row, local_db):

    if isinstance(local_db, str):
        local_table = get_local_packages(local_db)
    elif isinstance(local_db, Table):
        local_table = local_db
    else:
        raise ValueError("local_db must be either Table or path to DB file")

    if type(new_row) != Row:
        raise ValueError("pkg_entry must be an astropy.table.Row object")

    if new_row["name"] in local_table["name"]:
        ii = npwhere(local_table["name"] == new_row["name"])[0][0]
        fmt = '%Y-%m-%d'
        local_date = dt.datetime.strptime(local_table[ii]["date_modified"], fmt)
        new_date   = dt.datetime.strptime(new_row["date_modified"], fmt)

        if new_date >= local_date:

            dic = {col : local_table[ii][col] for col in local_table.colnames}
            dic["name"] = dic["name"] + "_" + dic["date_modified"]
            new_tbl = Table(names=[col for col in dic],
                            data=[[dic[col]] for col in dic])
            new_tbl_2 = Table(names=[col for col in new_row.colnames],
                              data=[[new_row[col]] for col in new_row.colnames])
            tbl_to_add = vstack([new_tbl, new_tbl_2])

            local_table.remove_row(ii)
        else:
            dic = {col: new_row[col] for col in new_row.colnames}
            dic["name"] = dic["name"] + "_" + dic["date_modified"]
            tbl_to_add = Table(names=[col for col in dic],
                               data=[[dic[col]] for col in dic])

    else:
        tbl_to_add = Table(names=[col for col in new_row.colnames],
                           data=[[new_row[col]] for col in new_row.colnames])

    new_local_table = vstack([local_table, tbl_to_add])
    new_local_table.meta = local_table.meta

    return new_local_table


def change_table_entry(tbl, col_name, old_val, new_val):

    offending_col = list(tbl[col_name].data)

    for ii in npwhere(old_val in offending_col)[0]:
        offending_col[ii] = new_val

    fixed_col = Column(name=col_name, data=offending_col)

    ii = npwhere(nparray(tbl.colnames) == col_name)[0][0]
    tbl.remove_column(col_name)
    tbl.add_column(fixed_col, index=ii)

    return tbl
