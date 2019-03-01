# import os
# import pytest
#
# import simcado as sim
#
# from simcado.tests.mocks.py_objects.source_objects import _image_source
#
# MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
#                                          "../mocks/MICADO_SCAO_WIDE/"))
# sim.rc.__search_path__ += [MOCK_PATH]
#
#
# ################################################################################
# # Source mocks
#
#
# @pytest.fixture(scope="function")
# def image_source():
#     return _image_source()


#
# import os
# import inspect
# import pytest
#
# from synphot import SpectralElement
#
# from simcado import UserCommands
# from simcado.optics import imager2 as imager
# import simcado as sim
#
#
# def mock_dir():
#     cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
#     rel_dirname = "mocks/MICADO_SCAO_WIDE/"
#
#     return os.path.abspath(os.path.join(cur_dirname, rel_dirname))
#
#
# MOCK_DIR = mock_dir()
# sim.rc.__search_path__.insert(0, MOCK_DIR)
#
#
# @pytest.fixture(scope="class")
# def opt_empty():
#     cmd = UserCommands()
#     opt = imager.Imager(cmd)
#
#     return opt
#
#
# @pytest.fixture(scope="class")
# def opt_scao_wide():
#     fname = os.path.join(MOCK_DIR, "mock_MICADO_SCAO_WIDE.config")
#     cmd = UserCommands(sim_data_dir=MOCK_DIR, filename=fname)
#     opt = imager.Imager(cmd)
#
#     return opt
#
#
# @pytest.mark.usefixtures("opt_scao_wide")
# class TestImagerInit:
#     def test_initialises_with_nothing(self):
#         opt = imager.Imager()
#         assert isinstance(opt, imager.Imager)
#
#     def test_accepts_usercommands(self, opt_scao_wide):
#         assert opt_scao_wide.cmds["INST_FILTER_TC"] == "Ks"
#
#
# @pytest.mark.usefixtures("opt_scao_wide", "opt_empty")
# class TestImagerSurfacesAttr:
#     def test_returns_empty_table_when_opt_is_empty(self, opt_empty):
#         srf_table = opt_empty.surfaces
#         assert len(srf_table) == 0
#
#     def test_returns_full_table_for_existing_table_files(self, opt_scao_wide):
#         srf_table = opt_scao_wide.surfaces
#         assert len(srf_table) == 19
#
#
# class TestMakeSurfacesTable:
#     def test_returns_empty_table_for_no_filenames(self):
#         surf_tbl = imager.make_surfaces_table()
#         assert len(surf_tbl) == 0
#
#     def test_returns_single_table_when_only_one_filename_is_passed(self):
#         files = ["EC_mirrors_ELT.tbl"]
#         surf_tbl = imager.make_surfaces_table(files)
#         assert len(surf_tbl) == 5
#         assert "Temp" in surf_tbl.colnames
#
#     def test_returns_combined_table(self):
#         files = ["EC_mirrors_ELT.tbl",
#                  "EC_mirrors_SCAO_relay.tbl",
#                  "EC_mirrors_MICADO_Wide.tbl"]
#         surf_tbl = imager.make_surfaces_table(files)
#         assert len(surf_tbl) == 19
#         assert "Temp" in surf_tbl.colnames
#
#     def test_returns_none_for_bogus_table(self):
#         surf_tbl = imager.make_surfaces_table(["bogus.tbl"])
#         assert len(surf_tbl) == 0
#
#     def test_ignores_tables_which_dont_exist_but_doesnt_throw_error(self):
#         files = ["EC_mirrors_ELT.tbl",
#                  "bogus.tbl"]
#         surf_tbl = imager.make_surfaces_table(files)
#         assert len(surf_tbl) == 5
#         assert "Temp" in surf_tbl.colnames
#
#
# @pytest.mark.usefixtures("opt_scao_wide")
# class TestMakeSpectralCurveFromFile:
#     def test_throw_exception_when_file_doesnt_exist(self, opt_scao_wide):
#         with pytest.raises(ValueError):
#             imager.import_spectral_curve_from_file("bogus.dat")
#
#     def test_reads_ok_for_existing_file(self):
#         file = "TC_filter_Ks.dat"
#         curve = imager.import_spectral_curve_from_file(file)
#         assert type(curve) == SpectralElement
#
#     def test_reads_reflectivity_if_exists(self):
#         file = "TER_dichroic.dat"
#         curve = imager.import_spectral_curve_from_file(file,
#                                                        val_name="reflection")
#         assert type(curve) == SpectralElement
#
#     def test_raises_error_if_colname_doesnt_exist(self):
#         pass
#
#
#
