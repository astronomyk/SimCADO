import os
import simcado as sim


class TestBasicLoading:

    def test_pkg_dir_in_search_path(self):
        assert sim.__pkg_dir__ in sim.__search_path__

    def test_data_dir_in_search_path(self):
        assert sim.__data_dir__ in sim.__search_path__

    def test_data_dir_in_pkg_dir(self):
        assert os.path.commonpath([sim.__pkg_dir__,
                                   sim.__data_dir__]) == sim.__pkg_dir__

    def test_rc_file_is_read_in(self):
        assert "SIM_LOGGING" in sim.__rc__

    def test_has_version_info(self):
        assert sim.__version__
