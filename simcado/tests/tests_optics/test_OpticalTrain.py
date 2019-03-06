import os
import pytest

import numpy as np

import simcado as sim
from simcado.optics.fov_manager import FOVManager
from simcado.optics import fov_manager as fov_mgr
from simcado.optics.image_plane import ImagePlane
from simcado.optics.optical_train import OpticalTrain
from simcado.optics.optics_manager import OpticsManager
from simcado.utils import find_file
from simcado.commands.user_commands2 import UserCommands
from simcado.optics.effects import TERCurve

from simcado.tests.mocks.py_objects.effects_objects import _mvs_effects_list
from simcado.tests.mocks.py_objects.yaml_objects import \
    _usr_cmds_min_viable_scope
from simcado.tests.mocks.py_objects.source_objects import _image_source, \
    _table_source

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
PLOTS = True

FILES_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/files/"))
YAMLS_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          "../mocks/yamls/"))

sim.rc.__search_path__ += [FILES_PATH, YAMLS_PATH]


def _basic_cmds():
    return UserCommands(filename=find_file("CMD_mvs_cmds.config"))


@pytest.fixture(scope="function")
def cmds():
    return _basic_cmds()


@pytest.fixture(scope="function")
def tbl_src():
    return _table_source()


@pytest.fixture(scope="function")
def im_src():
    return _image_source()




@pytest.mark.usefixtures("cmds")
class TestInit:
    def test_initialises_with_nothing(self):
        assert isinstance(OpticalTrain(), OpticalTrain)

    def test_initialises_with_basic_commands(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt, OpticalTrain)

    def test_has_observation_dict_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert len(opt.observation_dict) != 0

    def test_has_optics_manager_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.optics_manager, OpticsManager)

    def test_has_fov_manager_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        print(opt.fov_manager)
        assert isinstance(opt.fov_manager, FOVManager)

    def test_has_image_plane_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert isinstance(opt.image_plane, ImagePlane)

    def test_has_yaml_dict_object_after_initialising(self, cmds):
        opt = OpticalTrain(cmds=cmds)
        assert len(opt.yaml_dicts) == 4


@pytest.mark.usefixtures("cmds", "im_src", "tbl_src")
class TestObserve:
    def test_observe_works_for_table(self, cmds, tbl_src):
        opt = OpticalTrain(cmds)
        opt.observe(tbl_src)

        if PLOTS:
            plt.imshow(opt.image_plane.image.T, origin="lower", norm=LogNorm())
            plt.show()

    def test_observe_works_for_image(self, cmds, im_src):
        opt = OpticalTrain(cmds)
        opt.observe(im_src)

        if PLOTS:
            plt.imshow(opt.image_plane.image.T, origin="lower", norm=LogNorm())
            plt.show()

    def test_observe_works_for_source_distributed_over_several_fovs(self, cmds,
                                                                    im_src):
        print(im_src.fields)
        cmds["SIM_DETECTOR_PIX_SCALE"] = 0.1
        opt = OpticalTrain(cmds)
        opt.observe(im_src)

        if PLOTS:
            plt.imshow(opt.image_plane.image.T, origin="lower", norm=LogNorm())
            plt.show()

