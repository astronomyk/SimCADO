"""
End-to-end simulator for MICADO on the E-ELT
============================================


"""

from .version import version as __version__

from .detector import *
from .source import *
from .optics import *
from .psf import *
from .spatial import *
from .commands import *
from .spectral import *
#from .SpectralGrating import *
from .simulation import *
from .utils import *
