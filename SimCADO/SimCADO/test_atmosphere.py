import unittest
import AtmosphereModel as am


class TestRefraction(unittest.TestCase):

    known_values = ((2.2, 75.7114262),
                    (1.4, 75.8414404))

    def test_refraction(self):
        """atmospheric_refraction should (almost) reproduce known values"""
        for lam, value in self.known_values:
            result = am.atmospheric_refraction(lam)
            self.assertAlmostEqual(result, value)


    def test_zenith(self):
        """atmospheric_refraction should return zero at the zenith"""

        result = am.atmospheric_refraction(2.2, z0=0)
        self.assertEqual(result, 0.)

    def test_altitude(self):
        """atmospheric_refraction should be larger at lower altitudes"""
        res_low = am.atmospheric_refraction(2.2, h=0)
        res_high = am.atmospheric_refraction(2.2, h=3064)

        self.assertGreater(res_low, res_high)

    def test_wavelength(self):
        """atmospheric_refraction should be smaller for longer wavelengths"""

        res_blue = am.atmospheric_refraction(1.2)
        res_red = am.atmospheric_refraction(2.2)

        self.assertLess(res_red, res_blue)

if __name__ == "__main__":
    unittest.main()
