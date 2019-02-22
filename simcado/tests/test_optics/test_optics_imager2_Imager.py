import os
import pytest

from astropy.io import fits

import simcado as sim
from simcado.optics import imager as opt
from simcado.optics.effects import *

from simcado.tests.mocks.py_objects.yaml_objects import \
    _atmo_yaml_dict, _detector_yaml_dict, _inst_yaml_dict

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
sim.rc.__search_path__ += [MOCK_PATH]


@pytest.fixture(scope="function")
def atmo_yaml_dict():
    return _atmo_yaml_dict()


@pytest.fixture(scope="function")
def inst_yaml_dict():
    return _inst_yaml_dict()


@pytest.fixture(scope="function")
def detector_yaml_dict():
    return _detector_yaml_dict()


@pytest.mark.usefixtures("atmo_yaml_dict")
class TestMakeEffect:
    def test_it_creates_an_effects_object(self, atmo_yaml_dict):
        effdic = atmo_yaml_dict["effects"][0]
        properties = atmo_yaml_dict["properties"]
        effect = opt.make_effect(effdic, **properties)

        assert isinstance(effect, GaussianDiffractionPSF)
        assert effect.meta["diameter"] == 39


@pytest.mark.usefixtures("atmo_yaml_dict")
class TestOpticalElementInit:
    def test_initialised_with_nothing(self):
        assert isinstance(opt.OpticalElement(), opt.OpticalElement)

    def test_initialised_with_yaml_dict(self, atmo_yaml_dict):
        opt_el = opt.OpticalElement(atmo_yaml_dict)
        assert isinstance(opt_el, opt.OpticalElement)
        assert isinstance(opt_el.effects[0], GaussianDiffractionPSF)


@pytest.mark.usefixtures("detector_yaml_dict")
class TestOpticalElementGetZOrderEffects:
    @pytest.mark.parametrize("z_orders, n", [(0, 2), (100, 1), ([200, 299], 1)])
    def test_returns_the_effects_with_z_values(self, z_orders, n,
                                               detector_yaml_dict):
        opt_el = opt.OpticalElement(detector_yaml_dict)
        assert len(opt_el.get_z_order_effects(z_orders)) == n


@pytest.mark.usefixtures("detector_yaml_dict")
class TestOpticalElementSurfaceListProperty:
    def test_returns_empty_list_if_no_surface_list_given(self):
        pass


@pytest.mark.usefixtures("detector_yaml_dict", "inst_yaml_dict")
class TestOpticsManager:
    def test_initialises_with_nothing(self):
        assert isinstance(opt.OpticsManager(), opt.OpticsManager)

    def test_initialises_yaml_dict(self, detector_yaml_dict):
        opt_man = opt.OpticsManager(detector_yaml_dict)
        assert isinstance(opt_man, opt.OpticsManager)

    def test_initialises_yaml_dict(self, detector_yaml_dict, inst_yaml_dict):
        opt_man = opt.OpticsManager([detector_yaml_dict, inst_yaml_dict])
        print(opt_man)
        assert isinstance(opt_man, opt.OpticsManager)

    def test_has_effects_loaded(self, detector_yaml_dict):
        opt_man = opt.OpticsManager([detector_yaml_dict])
        # print(opt_man.optical_elements[1])
        assert isinstance(opt_man.optical_elements[1], opt.OpticalElement)
        assert isinstance(opt_man.optical_elements[1].effects[0], Effect)


@pytest.mark.usefixtures("detector_yaml_dict")
class TestOpticsManagerImagePlaneHeader:
    def test_makes_image_plane_header_correctly(self, detector_yaml_dict):
        opt_man = opt.OpticsManager(detector_yaml_dict)
        opt_man.meta["SIM_DETECTOR_PIX_SCALE"] = 0.004
        print(opt_man)
        assert isinstance(opt_man.image_plane_header, fits.Header)
































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
#
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
