import os
from copy import deepcopy
import pytest
from astropy.table import Table, Row
import telescopy.data.database as tdb
import telescopy.default_keywords as dkeys


def test_get_local_pkgs_returns_Table_if_local_DB_file_exists():
    local_db_path = os.path.join(dkeys.PKG_DIR,
                                 dkeys.INST_PKG_LOCAL_PATH,
                                 dkeys.INST_PKG_LOCAL_DB_NAME)
    assert type(tdb.get_local_packages(local_db_path)) == Table

def test_get_local_pkgs_throws_error_if_local_DB_file_path_is_bogus():
    local_db_path = "bogus.txt"
    with pytest.raises(ValueError):
        tdb.get_local_packages(local_db_path)


def test_get_server_packages_throws_error_on_wrong_path():
    svr_db_url = "www.my-server.bogus"
    with pytest.raises(ValueError):
        tdb.get_server_packages(svr_db_url)

def test_get_server_packages_returns_Table_if_path_correct():
    svr_path = os.path.join(dkeys.INST_PKG_SERVER_PATH,
                            dkeys.INST_PKG_SERVER_DB_NAME)
    assert type(tdb.get_server_packages(svr_path)) == Table


def test_rename_server_table_throws_exception_when_not_enough_column_names():
    svr_table = Table(data=[[1,1],[1,1]])
    svr_table.meta["comments"] = ["# name"]

    with pytest.raises(Exception):
        tdb.rename_table_colnames(svr_table)

def test_rename_server_table_renames_when_num_cols_equals_num_names():
    svr_table = Table(data=[[1,1],[1,1]])
    svr_table.meta["comments"] = ["# name other"]

    renamed_tbl = tdb.rename_table_colnames(svr_table)
    assert renamed_tbl.colnames[0] == "name"
    assert renamed_tbl.colnames[1] == "other"


def test_check_package_exists_returns_True_for_package_name():
    assert tdb.check_package_exists("test_package") == True

def test_check_package_exists_returns_exception_for_bogus_pkg_name():
    with pytest.raises(ValueError):
        tdb.check_package_exists("bogus")


def test_get_server_package_path_returns_url_with_right_data():
    svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                    ["test_package.zip"]])
    assert tdb.get_server_package_path("test_package", svr_table) == \
                                                            "test_package.zip"

def test_get_server_package_path_returns_exception_with_wrong_data():
    svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                    ["test_package.zip"]])
    with pytest.raises(ValueError):
        assert tdb.get_server_package_path("bogus", svr_table)


def test_get_package_table_entry_returns_url_with_right_data():
    svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                    ["test_package.zip"]])
    assert tdb.get_package_table_entry("test_package", svr_table)["path"] == \
                                                            "test_package.zip"

def test_get_package_table_entry_returns_exception_with_wrong_data():
    svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                    ["test_package.zip"]])
    with pytest.raises(ValueError):
        assert tdb.get_package_table_entry("bogus", svr_table)

def test_get_package_table_entry_returns_newest_with_multiple_entries():
    pass

def test_get_package_table_entry_returns_oldest_with_multiple_entries():
    pass


def test_download_package_raise_error_when_pkg_not_in_DB():
    with pytest.raises(ValueError):
        tdb.download_package("bogus")

def test_download_package_raise_error_when_pkg_file_doesnt_exist():
    with pytest.raises(ValueError):
        tdb.download_package("non_existent_pkg")

# avoid a test that is dependent on the network
# ::todo add this to the integration test suite
# def test_download_package_no_error_when_package_exists_and_is_in_DB():
#     pass

def test_add_pkg_to_local_db_adds_row():
    local_table = Table(names=["name", "author", "date_added",
                               "date_modified", "path"],
                        data=[["test_package"], ["Kieran Leschinski"],
                              ["2018-11-09"], ["2018-11-09"],
                              ["test_package.zip"]])
    len_local_table = len(local_table)
    pkg_entry = local_table[0]
    new_local_table = tdb.add_pkg_to_local_db(pkg_entry, local_table)

    assert len(new_local_table) == len_local_table + 1

def test_add_pkg_to_load_db_throws_exception_if_pkg_entry_isnt_Row():
    local_table = Table(names=["name", "path"], data=[["test_package"],
                                                      ["test_package.zip"]])
    with pytest.raises(ValueError):
        tdb.add_pkg_to_local_db("hello world!", local_table)

def test_add_pkg_to_local_db_if_two_equal_names_older_is_renamed():
    local_table = Table(names=["name", "date_modified"],
                        data=[["test_package"], ["2018-11-09"]])
    pkg_entry = deepcopy(local_table[0])
    pkg_entry["date_modified"] = "2018-11-14"

    new_local_table = tdb.add_pkg_to_local_db(pkg_entry, local_table)

    assert new_local_table["name"][0] != "test_package"


def test_change_table_entry_string_changed_successfully():
    tbl = Table(names=["id", "name"], data=[[0], ["seb skelly"]])
    tbl = tdb.change_table_entry(tbl, "name", "seb skelly", "seb skelly rocks")
    assert tbl[0]["name"] == "seb skelly rocks"

#def test_change_table_en




