import numpy as np

from .data_container import DataContainer


class LinearityCurve(DataContainer):
    def __init__(self, **kwargs):
        super(LinearityCurve, self).__init__(**kwargs)
        self.validate("LINEARIT")

    def apply_to(self, image, **kwargs):
        if not isinstance(image, np.ndarray):
            raise ValueError("`image` accepts only 2D-arrays ")

        out_array = np.interp(image.flatten(),
                              self.table["real_flux"],
                              self.table["detected_flux"]).reshape(image.shape)
        out_array = out_array.astype(np.float32)

        return out_array


class TERCurve(DataContainer):
    def __init__(self, **kwargs):
        super(TERCurve, self).__init__(**kwargs)
        self.validate("TERCURVE")

    def apply_to(self, source_spectrum, **kwargs):
        return source_spectrum


class ConstantPSF(DataContainer):
    def __init__(self, **kwargs):
        super(TERCurve, self).__init__(**kwargs)
        self.validate("CONSTPSF")

    def apply_to(self, source, **kwargs):
        pass












