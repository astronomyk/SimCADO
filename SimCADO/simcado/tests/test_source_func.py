import numpy as np
import simcado as sim

def test_BV_to_spec_type():
    """
    Test the border cases:
    - O1V for B-V < -0.5,
    - M9V for B-V > 3
    Test for string and list

    """
    blue       = sim.source.BV_to_spec_type(-0.5)
    red        = sim.source.BV_to_spec_type(3)
    list_BV    = sim.source.BV_to_spec_type([0, 1, 2])
    single_BV  = sim.source.BV_to_spec_type(1.0)

    assert blue == "O0V"

    assert red == "M9V"

    assert type(list_BV) == list
    assert list_BV == ['A4V', 'K1V', 'M6V']

    assert single_BV == "K1V"


def test_mag_to_photons():
    """
    Test
    - V=0,
    - V=20,
    - Ks=0,
    - Ks=30
    """
    v0  = mag_to_photons("V", 0)
    v20 = mag_to_photons("V", 20)
    k0  = mag_to_photons("V", 0)
    k30 = mag_to_photons("V", 30)

    # got from cfa website
    assert abs((v0 - 8786488925.436462) / v0) < 0.05
