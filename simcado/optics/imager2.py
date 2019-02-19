import os
import warnings
from copy import deepcopy

from astropy.io import ascii as ioascii
from astropy.table import Table, vstack
from astropy import units as u

from synphot import SpectralElement
from synphot.models import Empirical1D

from ..utils import find_file


class OpticalTrain():
    def load_optical_elements(self, filename):

        self.optical_elements = []


    def make_radiometry_table(self, filename):

        self.radiometry_table = Table



class Imager:
    pass