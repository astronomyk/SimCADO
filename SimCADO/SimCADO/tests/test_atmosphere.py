import unittest
import pytest
import numpy as np
from ..AtmosphereModel import *

class TestRefraction(unittest.TestCase):
    """This class uses unittest"""

    known_values = ((2.2, 75.7114262),
                    (1.4, 75.8414404))

    def test_refraction(self):
        """atmospheric_refraction should (almost) reproduce known values"""
        for lam, value in self.known_values:
            result = atmospheric_refraction(lam)
            self.assertAlmostEqual(result, value)


    def test_zenith(self):
        """atmospheric_refraction should return zero at the zenith"""

        result = atmospheric_refraction(2.2, z0=0)
        self.assertEqual(result, 0.)

    def test_altitude(self):
        """atmospheric_refraction should be larger at lower altitudes"""
        res_low = atmospheric_refraction(2.2, h=0)
        res_high = atmospheric_refraction(2.2, h=3064)

        self.assertGreater(res_low, res_high)

    def test_wavelength(self):
        """atmospheric_refraction should be smaller for longer wavelengths"""

        res_blue = atmospheric_refraction(1.2)
        res_red = atmospheric_refraction(2.2)

        self.assertLess(res_red, res_blue)

class TestPytest:
    """This class uses py.test"""
    known_values = ((2.2, 75.7114262),
                    (1.4, 75.8414404))

    def test_refraction(self):
        """atmospheric_refraction should (almost) reproduce known values"""
        for lam, value in self.known_values:
            result = atmospheric_refraction(lam)
            assert np.isclose(result, value)

    def test_zenith(self):
        """atmospheric_refraction should return zero at the zenith"""

        result = atmospheric_refraction(2.2, z0=0)
        assert result == 0.

    def test_altitude(self):
        """atmospheric_refraction should be larger at lower altitudes"""
        res_low = atmospheric_refraction(2.2, h=0)
        res_high = atmospheric_refraction(2.2, h=3064)

        assert res_low > res_high

    def test_wavelength(self):
        """atmospheric_refraction should be smaller for longer wavelengths"""

        res_blue = atmospheric_refraction(1.2)
        res_red = atmospheric_refraction(2.2)

        assert res_red < res_blue

if __name__ == "__main__":
    # unittest.main()
    pytest.main("-v")
