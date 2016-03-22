"""
End-to-end simulator for MICADO on the E-ELT
============================================


"""

from .version import version as __version__

from .Detector import *
from .LightObject import *
from .OpticalTrain import *
from .AtmosphereModel import *
from .PSFCube import *
from .PlaneEffect import *
from .UserCommands import *
from .SpectralCurve import *
from .SpectralGrating import *
from .SimulationRun import *
from .utils import *
