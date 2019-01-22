# 1 read in the tables
# 2 read in curves from the set of unique files
# 3 create a dictionary of curves
#
import pytest
import os
import inspect

import numpy as np
from astropy.table import Table
from astropy.io import ascii as ioascii

from synphot import SpectralElement, SourceSpectrum

from simcado.optics import radiometry as opt_rad
from simcado.optics import surface as opt_surf
from simcado import utils
import simcado as sim


def mock_dir():
    cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
    rel_dirname = "mocks/MICADO_SCAO_WIDE/"
    return os.path.abspath(os.path.join(cur_dirname, rel_dirname))


MOCK_DIR = mock_dir()
sim.rc.__search_path__.insert(0, MOCK_DIR)


@pytest.fixture(scope="module")
def input_tables():
    filenames = ["EC_mirrors_ELT.tbl",
                 "EC_mirrors_SCAO_relay.tbl",
                 "EC_mirrors_MICADO_Wide.tbl"]

    return [os.path.join(MOCK_DIR, fname) for fname in filenames]


@pytest.mark.usefixtures("input_tables")
class TestInit:
    def test_initialises_with_no_input(self):
        rt = opt_rad.RadiometryTable()
        assert isinstance(rt, opt_rad.RadiometryTable)
        assert rt.table is None

    def test_initialises_with_single_table(self, input_tables):
        rt = opt_rad.RadiometryTable([input_tables[0]])
        assert isinstance(rt, opt_rad.RadiometryTable)
        assert len(rt.table) == 5

    def test_initialises_with_list_of_tables(self, input_tables):
        rt = opt_rad.RadiometryTable(input_tables)
        assert isinstance(rt, opt_rad.RadiometryTable)
        assert len(rt.table) == 19


@pytest.mark.usefixtures("input_tables")
class TestAddSurfaceList:
    def test_append_single_table_from_filename(self, input_tables):
        rad_table = opt_rad.RadiometryTable()
        rad_table.add_surface_list(input_tables[0])
        assert len(rad_table.table) == 5

    def test_combine_two_tables_from_filename(self, input_tables):
        rad_table = opt_rad.RadiometryTable()
        rad_table.add_surface_list([input_tables[0], input_tables[1]])
        assert len(rad_table.table) == 6

    def test_combine_list_of_filename(self, input_tables):
        rad_table = opt_rad.RadiometryTable()
        rad_table.add_surface_list(input_tables)
        assert len(rad_table.table) == 19


@pytest.mark.usefixtures("input_tables")
class TestCombineTables:
    def test_adds_two_tables(self):
        tblA = Table(names=["colA", "colB"], data=[[0, 1], [0, 1]])
        tblB = Table(names=["colA", "colB"], data=[[2, 3], [2, 3]])
        tblC = opt_rad.combine_tables(tblB, tblA)
        assert np.all(tblC["colB"] == np.array([0, 1, 2, 3]))

    def test_adds_single_table(self):
        tblA = Table(names=["colA", "colB"], data=[[0, 1], [0, 1]])
        tblC = opt_rad.combine_tables(tblA)
        assert np.all(tblC["colA"] == np.array([0, 1]))

    def test_adds_three_tables_to_old_table(self):
        tblA = Table(names=["colA", "colB"], data=[[0, 1], [0, 1]])
        tblB = Table(names=["colA", "colB"], data=[[2, 3], [2, 3]])
        tblC = Table(names=["colA", "colB"], data=[[4, 5], [4, 5]])
        tblD = Table(names=["colA", "colB"], data=[[6, 7], [6, 7]])
        tblE = opt_rad.combine_tables([tblB, tblC, tblD], tblA)
        assert np.all(tblE["colA"] == np.arange(8))

    def test_adds_table_from_filename_to_nothing(self, input_tables):
        tblA = ioascii.read(input_tables[0])
        tblC = opt_rad.combine_tables(tblA)
        assert len(tblC) == 5

    def test_adds_table_from_filename_to_table_object(self, input_tables):
        tblA = ioascii.read(input_tables[0])
        tblB = input_tables[1]
        tblC = opt_rad.combine_tables(tblB, tblA)
        assert len(tblC) == 6

    def test_adds_table_from_filename_to_table_from_file(self, input_tables):
        tblA = input_tables[0]
        tblB = input_tables[1]
        tblC = opt_rad.combine_tables(tblB, tblA)
        assert len(tblC) == 6

    def test_adds_3_tables_from_filename_to_nothing(self, input_tables):
        tblC = opt_rad.combine_tables(input_tables)
        assert len(tblC) == 19

    def test_prepend_table(self):
        tblA = Table(names=["colA", "colB"], data=[[0, 1], [0, 1]])
        tblB = Table(names=["colA", "colB"], data=[[2, 3], [2, 3]])
        tblC = opt_rad.combine_tables(tblB, tblA, prepend=True)
        assert np.all(tblC["colB"] == np.array([2, 3, 0, 1]))


@pytest.mark.usefixtures("input_tables")
class TestMakeSurfaceFromRow:
    def test_return_none_from_empty_row(self, input_tables):
        tblA = ioascii.read(input_tables[0])
        surf = opt_rad.make_surface_from_row(tblA[0])
        assert isinstance(surf, opt_surf.SpectralSurface)

    def test_surface_has_processed_ter_filename_in_row(self, input_tables):
        tblA = ioascii.read(input_tables[0])
        surf = opt_rad.make_surface_from_row(tblA[0])
        assert isinstance(surf.transmission, SpectralElement)
        assert isinstance(surf.reflection, SpectralElement)
        assert isinstance(surf.emissivity, SpectralElement)
        assert isinstance(surf.emission, SourceSpectrum)


class TestGetRealColname:
    @pytest.mark.parametrize("name, colnames", [
                             ("yahoo", ["Yahoo", "Bogus"]),
                             ("yahoo", ["yahoo", "Bogus"]),
                             ("yahoo", ["YAHOO", "Bogus"])])
    def test_returns_real_name(self, name, colnames):
        assert opt_rad.real_colname(name, colnames) == colnames[0]


@pytest.mark.usefixtures("input_tables")
class TestMakeSurfaceDictFromTable:
    def test_return_dict_from_table(self, input_tables):
        tbl = ioascii.read(input_tables[0])
        surf_dict = opt_rad.make_surface_dict_from_table(tbl)
        assert isinstance(surf_dict, dict)


