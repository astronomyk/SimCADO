import unittest
import pytest
import numpy as np
from ..SpectralCurve import *

class TestNormalize(unittest.TestCase):
    """This class uses unittest"""

    def test_normalize_integral(self):
        """Test normalization of unit transmission"""
        lam = np.linspace(1.0, 2.0, 101)
        val = 2.
        normval = 1.

        ucurve = UnityCurve(lam=lam, val=val)
        ucurve.normalize(val=normval, mode='integral')
        unitycheckval = np.sum(ucurve.val[1:] * \
                               (ucurve.lam[1:] - ucurve.lam[:-1]))

        self.assertAlmostEqual(unitycheckval, normval)

    def test_normalize_maximum(self):
        """Test maximum normalization of unit transmission"""
        lam = np.linspace(1.0, 2.0, 101)
        val = 2.
        normval = 1.

        unitycurve = UnityCurve(lam=lam, val=val)
        unitycurve.normalize(val=normval, mode='maximum')
        unitycheckval = np.max(unitycurve.val)
        self.assertAlmostEqual(unitycheckval, normval)

    def test_normalize_error(self):
        """normalize should raise ValueError with wrong mode"""
        lam = np.linspace(1.0, 2.0, 101)
        unitycurve = UnityCurve(lam=lam, val=2.)
        self.assertRaises(ValueError,
                          unitycurve.normalize, mode="Wrong mode")

if __name__ == "__main__":
    #unittest.main()
    pytest.main("-v")
