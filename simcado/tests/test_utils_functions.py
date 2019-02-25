"""Unit tests for module simcado.utils"""

import pytest
import numpy as np
from astropy import wcs
from astropy.io import ascii as ioascii, fits

from simcado.utils import parallactic_angle, deriv_polynomial2d
from simcado.utils import find_file, has_needed_keywords
from simcado.utils import airmass2zendist, zendist2airmass
from simcado.utils import convert_table_comments_to_dict

from simcado import rc


class TestFindFile:
    """Tests of function simcado.utils.find_file"""

    def test_fails_if_filename_not_a_string(self):
        # python 3.6: TypeError
        # python 3.4, 3.5: AttributeError (change in os.path.isabs)
        with pytest.raises((TypeError, AttributeError)):
            find_file(1.2, rc.__search_path__)

    def test_passes_if_file_exists(self):
        filename = 'utils.py'
        assert find_file(filename, rc.__search_path__)

    def test_fails_if_file_doesnt_exist(self):
        filename = 'utils987654.pz'
        assert find_file(filename, rc.__search_path__) is None

    def test_ignores_none_objects_in_search_path_list(self):
        filename = 'utils.py'
        new_filename = find_file(filename, [None] + rc.__search_path__)
        assert filename in new_filename


class TestAirmassZendist:
    """Tests conversion between airmass and zenith distance"""

    def test_airmass2zendist_pass_for_known_quanities_AM_1_equals_ZD_0(self):
        assert np.allclose(airmass2zendist(1.0), 0)

    def test_pass_for_known_quanities_AM_2_equals_ZD_sqrt2(self):
        assert np.allclose(airmass2zendist(np.sqrt(2)), 45)

    def test_zendist2airmass_pass_for_known_quanities_ZD_0_equals_AM_1(self):
        assert np.allclose(zendist2airmass(0), 1.0)

    def test_zendist2airmass_pass_for_known_quanities_ZD_60_equals_ZD_2(self):
        assert np.allclose(zendist2airmass(60), 2.0)

    def test_zendist2airmass_undoes_exactly_what_airmass2zendist_does(self):
        airmass = 1.78974234
        assert np.allclose(zendist2airmass(airmass2zendist(airmass)),
                           airmass)

    def test_airmass2zendist_undoes_exactly_what_zendist2airmass_does(self):
        zendist = 12.31334
        assert np.allclose(airmass2zendist(zendist2airmass(zendist)),
                           zendist)


class TestParallacticAngle:
    """Tests of function simcado.utils.parallactic_angle"""

    def test_parallactic_angle_negative_east_of_meridian(self):
        assert parallactic_angle(-1, 0, -24) < 0

    def test_parallactic_angle_positive_west_of_meridian(self):
        assert parallactic_angle(1, 0, -24) > 0

    def test_parallactic_angle_zero_on_meridian(self):
        assert parallactic_angle(0, 0, 24) == 0

    def test_specific_example_from_Ball_1908(self):
        """Test: Example from Ball (1908), p.92"""
        ha = -3.                 # 3 hours east
        de = 38 + 9/60.          # decl 38d09m
        lat = 53 + 23/60.        # lat  53d23m
        eta0 = - (48 + 41/60.)   # result -48d41m

        eta = parallactic_angle(ha, de, lat)

        # should agree to within 1 arcmin
        assert np.allclose(eta, eta0, atol=1/60.)

    def test_setting_object_on_the_equator_is_90_minus_latitude(self):
        """
        For a setting object on the equator, the parallactic angle
        is 90 - lat
        """
        lat = np.random.rand(10) * 180 - 90
        pa = parallactic_angle(6, 0, lat)

        assert np.allclose(pa, 90. - lat)


class TestDerivPolynomial2D:
    """Tests of simcado.utils.deriv_polynomial2d"""

    def test_derivative_of_2D_polynomial_equal_to_analytical_derivative(self):
        from astropy.modeling.models import Polynomial2D

        ximg, yimg = np.meshgrid(np.linspace(-1, 1, 101),
                                 np.linspace(-1, 1, 101))
        poly = Polynomial2D(2, c0_0=1, c1_0=2, c2_0=3,
                            c0_1=-1.5, c0_2=0.4, c1_1=-2)
        # Expected values
        y_x = 2 + 6 * ximg - 2 * yimg
        y_y = -1.5 + 0.8 * yimg - 2 * ximg

        dpoly_x, dpoly_y = deriv_polynomial2d(poly)
        # Computed values
        y_x_test = dpoly_x(ximg, yimg)
        y_y_test = dpoly_y(ximg, yimg)

        assert np.allclose(y_x, y_x_test)
        assert np.allclose(y_y, y_y_test)


class TestConvertCommentsToDict:
    def test_converts_list_of_strings_to_dict_if_comments_in_table_meta(self):
        tbl = ioascii.read("""
                         # key1 : val 1 
                         # key2 : extra long entry
                         col1    col2
                         0       1 """)
        dic = convert_table_comments_to_dict(tbl)
        assert dic["key1"] == "val 1"
        assert len(dic) == 2

    def test_returns_none_if_comments_not_in_table_meta(self):
        tbl = ioascii.read("""col1    col2
                              0       1 """)
        dic = convert_table_comments_to_dict(tbl)
        assert dic is None

    def test_returns_input_if_conversion_doesnt_work(self):
        tbl_str = """
        # key1 : val 1 
        # 
        # key2
        col1    col2
        0       1 """
        tbl = ioascii.read(tbl_str)
        dic = convert_table_comments_to_dict(tbl)
        assert dic == tbl.meta["comments"]


class TestHasWcsKeys:
    def test_fails_if_header_does_not_have_all_keys(self):
        assert not has_needed_keywords(fits.Header())

    def test_passes_if_header_does_have_all_keys(self):
        hdr = wcs.WCS().to_header()
        hdr["NAXIS1"] = 100
        assert has_needed_keywords(hdr)

    def test_passes_if_header_does_have_all_keys_and_suffix(self):
        hdr = wcs.WCS(key="D").to_header()
        hdr["NAXIS1"] = 100
        assert has_needed_keywords(hdr, "D")

