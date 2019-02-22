import os
import pytest
from pytest import approx

from astropy import units as u

import simcado as sim
from simcado.optics.effects import DetectorList

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
sim.rc.__search_path__ += [MOCK_PATH]


class TestDetectorInit:
    def test_initialises_with_nothing(self):
        assert isinstance(DetectorList(), DetectorList)

    def test_initialises_with_filename(self):
        det_list = DetectorList(filename="FPA_array_layout.dat")
        assert isinstance(det_list, DetectorList)


class TestDetectorImagePlaneHeader:
    def test_header_is_sensical(self):
        det_list = DetectorList(filename="FPA_array_layout.dat")

        hdr_small = det_list.image_plane_header(1.5*u.mas)
        hdr_big = det_list.image_plane_header(4*u.mas)

        assert hdr_small["NAXIS1"] == approx(hdr_big["NAXIS1"], abs=2)
