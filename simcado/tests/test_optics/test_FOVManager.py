import os
import pytest

from astropy.io import fits

import simcado as sim
from simcado.optics import optical_train as opt
from simcado.optics.effects import *
from simcado.optics.fov_manager import FOVManager

from simcado.tests.mocks.py_objects.source_objects import _image_source
from simcado.tests.mocks.py_objects.yaml_objects import \
    _atmo_yaml_dict, _detector_yaml_dict, _inst_yaml_dict
from simcado.tests.mocks.py_objects.effects_objects import \
    _surf_list, _surf_list_empty, _filter_surface


class TestInit():
    def test_initialises_with_nothing(self):
        assert isinstance(FOVManager(), FOVManager)

    def test_initialises_with_list_of_effects(self):
        pass
