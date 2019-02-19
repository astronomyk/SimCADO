import numpy as np
from simcado.optics.effects.effects import Effect


class TestEffectInit:
    def test_initialises_with_no_input(self):
        assert isinstance(Effect(), Effect)

    def test_initalising_with_arrays_creates_table(self):
        eff = Effect(array_dict={"x": [-1, 0, 1], "y": [1, 0, 1],
                                 "flux": [1, 2, 3]})
        assert isinstance(eff, Effect)
        assert np.sum(eff.table["flux"]) == 6

    def test_has_method_apply_to(self):
        assert hasattr(Effect(), "apply_to")

    def test_has_method_fov_grid(self):
        assert hasattr(Effect(), "fov_grid")
