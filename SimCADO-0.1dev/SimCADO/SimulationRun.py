###############################################################################
# SimulationRun
#
# Deliverables
#
#
# Need 2 different UserCommands dictionaries as input
#  - Observation parameters (Alt, Az, Exptime, etc)
#  - Observatory parameters (area, instrument configuration, psfs, etc)
#
#
#
# Classes:
#
#
# Methods:
#
#
# 1. Read in UserCommands
#    Create UserCommands objects for: optical train and observing run.
#    The observing run will be used for science, atmo and mirror photons.
#
# 2. Create OpticalTrain, 
#    which in turn creates all the individual objects
#
# 3. Create LightObjects.
#    This in turn means creating a LightObject for every source of background
#    photons, e.g. for the atmosphere and every mirror.
#    If the cube is too large (>1GB), the LightObject for the science case will
#    need to be split into x lists (based on integrated pixel brightness) and
#    the simulation run x times with the same OpticalTrain.
#
# 4. To all the LightObjects do the following:
#    - Apply spectral effects (Transmission curves)
#    - Apply plane effects (distortion, rotation), where applicable
#    - Apply 3D PSF cube effects (Telescope PSF, ADC PSF), where applicable
#    - Collapse the cube into a 2D image
#    - Apply 2D PSF effects
#
# 5. Convert image to electrons and add the detector effects



######## Taken out of my IPython Notebook - not to be take seriously

import numpy as np
from astropy.io import fits
import sys

import PSFCube, SpectralCurve, LightObject2



