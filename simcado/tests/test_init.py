import sys
import os
import yaml
import simcado as sim
from .. import __pkg_dir__

class TestBasicLoading:

    def test_pkg_dir_in_search_path(self):
        assert sim.__pkg_dir__ in sim.__search_path__

    def test_data_dir_in_search_path(self):
        assert sim.__data_dir__ in sim.__search_path__

    def test_data_dir_in_pkg_dir(self):
        if sys.version_info >= (3,5):
            cpath = os.path.commonpath([sim.__pkg_dir__, sim.__data_dir__])
            assert cpath == sim.__pkg_dir__
        else:
            cpath = os.path.commonprefix([sim.__pkg_dir__, sim.__data_dir__])
            assert cpath == sim.__pkg_dir__

    def test_rc_file_is_read_in(self):
        assert "SIM_LOGGING" in sim.__rc__

    def test_has_version_info(self):
        assert sim.__version__


class TestRcFile:

    def test_rcfile_exists(self):
        assert os.path.exists(os.path.join(__pkg_dir__, ".simcadorc"))

    def test_rc_file_readable_by_yaml(self):
        with open(os.path.join(__pkg_dir__, ".simcadorc")) as rc_file:
            rc_dict = yaml.load(rc_file)
        assert isinstance(rc_dict, dict)
        assert len(rc_dict) > 0
