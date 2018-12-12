import pytest
from astropy.table import Table
import simcado.server.database as sim_db

_parent_path = "./test_downloads_dir/"
sim_db.set_up_local_package_directory(_parent_path, True)

# for path in local_dbs:
#     os.remove(path)
# os.removedirs(parent_path)


class TestGetLocalPackages:
    def test_returns_table_if_local_db_file_exists(self):
        local_dbs = sim_db.set_local_path_names(_parent_path)
        for path in local_dbs:
            db_tbl = sim_db.get_local_packages(path)
            assert isinstance(db_tbl, Table)

    def test_throws_error_if_local_db_file_path_is_bogus(self):
        local_db_path = "bogus.txt"
        with pytest.raises(ValueError):
            sim_db.get_local_packages(local_db_path)

    def test_empty_table_has_column_type_string(self):
        pass

class TestGetServerPackages:
    def test_throws_error_on_wrong_path(self):
        svr_db_url = "www.my-server.bogus"
        with pytest.raises(ValueError):
            sim_db.get_server_packages(svr_db_url)

    def test_returns_table_if_path_correct(self):
        svr_path = sim_db.SVR_INST_DB
        assert type(sim_db.get_server_packages(svr_path)) == Table
        svr_path = sim_db.SVR_PSF_DB
        assert type(sim_db.get_server_packages(svr_path)) == Table
        svr_path = sim_db.SVR_SRC_DB
        assert type(sim_db.get_server_packages(svr_path)) == Table


class TestRenameServerTable:
    def test_throws_exception_when_not_enough_column_names(self):
        svr_table = Table(data=[[1, 1], [1, 1]])
        svr_table.meta["comments"] = ["# name"]

        with pytest.raises(Exception):
            sim_db.rename_table_colnames(svr_table)

    def test_renames_when_num_cols_equals_num_names(self):
        svr_table = Table(data=[[1, 1], [1, 1]])
        svr_table.meta["comments"] = ["# name other"]

        renamed_tbl = sim_db.rename_table_colnames(svr_table)
        assert renamed_tbl.colnames[0] == "name"
        assert renamed_tbl.colnames[1] == "other"


class TestCheckPackageExists:
    def test_returns_true_for_package_name(self):
        assert sim_db.check_package_exists("test_package") is True

    def test_returns_exception_for_bogus_pkg_name(self):
        with pytest.raises(ValueError):
            sim_db.check_package_exists("bogus")


class TestGetServerPackagePath:
    def test_return_url_for_existing_package(self):
        svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                        ["test_package.zip"]])
        path = sim_db.get_server_package_path("test_package", svr_table)
        assert path == "test_package.zip"

    def test_returns_none_if_package_not_in_dbs(self):
        svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                        ["test_package.zip"]])
        path = sim_db.get_server_package_path("bogus", svr_table)
        assert path is None


class TestGetPackageTableEntry:
    def test_returns_url_with_right_data(self):
        svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                        ["test_package.zip"]])
        return_tbl = sim_db.get_package_table_entry("test_package", svr_table)
        assert return_tbl["path"] == "test_package.zip"

    def test_returns_exception_with_wrong_data(self):
        svr_table = Table(names=["name", "path"], data=[["test_package"],
                                                        ["test_package.zip"]])
        with pytest.raises(ValueError):
            assert sim_db.get_package_table_entry("bogus", svr_table)

    def test_returns_newest_with_multiple_entries(self):
        pass

    def test_returns_oldest_with_multiple_entries(self):
        pass


class TestDownloadPackage:
    def test_raise_error_when_pkg_not_in_db(self):
        with pytest.raises(ValueError):
            sim_db.download_package("bogus")

    def test_raise_error_when_pkg_file_doesnt_exist(self):
        with pytest.raises(ValueError):
            sim_db.download_package("non_existent_pkg")

# avoid a test that is dependent on the network
# ::todo add this to the integration test suite
# def test_download_package_no_error_when_package_exists_and_is_in_DB():
#     pass


class TestAddPackageToLocalDb:
    def test_adds_row(self):
        local_table = Table(names=["name", "author", "date_added",
                                   "date_modified", "path"],
                            data=[["test_package"], ["Kieran Leschinski"],
                                  ["2018-11-09"], ["2018-11-09"],
                                  ["test_package.zip"]])
        len_local_table = len(local_table)
        pkg_entry = local_table[0]
        new_local_table = sim_db.add_pkg_to_local_db(pkg_entry, local_table)

        assert len(new_local_table) == len_local_table + 1

    def test_throws_exception_if_pkg_entry_isnt_astropy_row_class(self):
        local_table = Table(names=["name", "path"], data=[["test_package"],
                                                          ["test_package.zip"]])
        with pytest.raises(ValueError):
            sim_db.add_pkg_to_local_db("hello world!", local_table)

    def test_existing_is_renamed_if_new_package_is_newer(self):
        local_table = Table(names=["name", "date_modified"],
                            data=[["test_package"], ["2018-11-09"]])
        pkg_entry = Table(names=["name", "date_modified"],
                          data=[["test_package"], ["2018-11-14"]])
        new_local_table = sim_db.add_pkg_to_local_db(pkg_entry[0], local_table)
        assert new_local_table["name"][0] == "test_package_2018-11-09"
        assert new_local_table["name"][1] == "test_package"

    def test_new_package_is_renamed_if_new_package_is_older(self):
        local_table = Table(names=["name", "date_modified"],
                            data=[["test_package"], ["2018-11-14"]])
        pkg_entry = Table(names=["name", "date_modified"],
                          data=[["test_package"], ["2018-11-09"]])
        new_local_table = sim_db.add_pkg_to_local_db(pkg_entry[0], local_table)
        assert new_local_table["name"][0] == "test_package"
        assert new_local_table["name"][1] == "test_package_2018-11-09"


class TestChangeTableEntry:
    def test_string_changed_successfully(self):
        tbl = Table(names=["id", "name"], data=[[0], ["seb skelly"]])
        tbl = sim_db.change_table_entry(tbl, "name", "seb skelly",
                                        "seb skelly rocks")
        assert tbl[0]["name"] == "seb skelly rocks"
