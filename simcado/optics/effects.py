import numpy as np

from .data_container import DataContainer

class Effect(DataContainer):

    def __init__(self):
        self.meta
        self.action_function
        self.children

    def apply_to(self, fov):


        return fov





class ConstantPSF(DataContainer):
    def __init__(self, **kwargs):
        super(ConstantPSF, self).__init__(**kwargs)
        self.validate("CONSTPSF")

    def apply_to(self, source, **kwargs):
        pass












