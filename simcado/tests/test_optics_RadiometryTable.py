# 1 read in the tables
# 2 read in curves from the set of unique files
# 3 create a dictionary of curves
#
import pytest
import os
import inspect


def mock_dir():
    cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
    rel_dirname = "mocks/MICADO_SCAO_WIDE/"

    return os.path.abspath(os.path.join(cur_dirname, rel_dirname))


MOCK_DIR = mock_dir()


@pytest.fixture(scope="module")
def input_tables():
    filenames = ["EC_mirrors_ELT.tbl",
                 "EC_mirrors_MICADO_Wide.tbl",
                 "EC_mirrors_SCAO_relay.tbl"]

    return [os.path.join(MOCK_DIR, fname) for fname in filenames]
