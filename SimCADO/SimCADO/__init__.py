"""
End-to-end simulator for MICADO on the E-ELT
============================================
"""
# Import all the modules to go under simcado.detector
from . import utils

from . import spectral
from . import spatial
from . import psf

from . import detector
from . import optics  
from . import commands
from . import source  

from . import optics_utils
from . import defaults

from .version import version as __version__

# import specific Classes from the modules to be accessible in the global 
# namespace
from .detector  import Detector, Chip
from .source    import Source
from .optics    import OpticalTrain
from .commands  import UserCommands

# don't import these ones just yet
#from .SpectralGrating  import * 
from . import simulation

