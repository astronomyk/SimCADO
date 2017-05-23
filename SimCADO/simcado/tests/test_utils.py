'''Unit tests for module simcado.utils'''

from simcado.utils import parallactic_angle
import numpy as np

class TestParallacticAngle():
    '''Tests of function simcado.utils.parallactic_angle'''

    def test_01(self):
        '''Test: parallactic angle negative east of meridian'''
        assert parallactic_angle(-1, 0, -24) < 0

    def test_02(self):
        '''Test: parallactic angle positive west of meridian'''
        assert parallactic_angle(1, 0, -24) > 0

    def test_03(self):
        '''Test: parallactic angle zero on meridian'''
        assert parallactic_angle(0, 0, 24) == 0

    def test_04(self):
        '''Test: Example from Ball (1908), p.92'''
        ha = -3                 # 3 hours east
        de = 38 + 9/60          # decl 38d09m
        lat = 53 + 23/60        # lat  53d23m
        eta0 = - (48 + 41/60)   # result -48d41m

        eta = parallactic_angle(ha, de, lat)

        # should agree to within 1 arcmin
        assert np.allclose(eta, eta0, atol=1/60)

    def test_05(self):
        '''Test parallactic angle

        For a setting object on the equator, the parallactic angle is 90 - lat'''
        lat = np.random.rand(10) * 180 - 90
        pa = parallactic_angle(6, 0, lat)

        assert np.allclose(pa, 90. - lat)
