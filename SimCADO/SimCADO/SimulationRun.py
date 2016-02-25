###############################################################################
# SimulationRun
#
# DESCIPTION
#
# File Structure needed to run the simulation
#
#
#
# Need 2 different UserCommands dictionaries as input
#  - Observation parameters (Alt, Az, Exptime, etc)
#  - Observatory parameters (area, instrument configuration, psfs, etc)
#  -
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

import sys, os
import warnings

import numpy as np

from astropy.io import fits
import astropy.units as u

try:
    import SimCADO.PSFCube as psf
    import SimCADO.SpectralCurve as sc
    import SimCADO.LightObject as lo
    import SimCADO.UserCommands as uc
    import SimCADO.utils as utils
except:
    import PSFCube as psf
    import SpectralCurve as sc
    import LightObject as lo
    import UserCommands as uc
    import utils

class Simulation(object):
    """
    """

    def __init__(self, cmds=None, filename=None):
        """
        """
        self.set_user_commands(filename)
        if  cmds  is not None: self.cmds = cmds

        self.check_defaults()


    def set_user_commands(self, user_filename=None):
        """
        - get master.config
        - read in all the other configs and add them to cmds
        - read in the file with user overrides

        Optional keywords:
        user_filename: path to the user.config file
        """
        self.cmds = utils.read_config("../user_commands/master.config")
        fnames = [self.cmds[key] for key in self.cmds.keys()]
        for fname in fnames:
            if os.path.exists(fname):
                self.cmds.update(utils.read_config(fname))
            else:
                warnings.warn(fname+" doesn't exist")

        if user_filename is not None:
            self.cmds.update(utils.read_config(user_filename))


    def read_optical_train(self, filename):
        # if we should use a previous optical system, import is
        pass

    def save_optical_train(self, filename, clobber=False):
        # save the current optical system to a FITS file
        pass

    def make_optical_train(self):
        # if we should use a previous optical system, import it
        # else generate a new one
        pass

    def read_source(self, filename):
        pass

    def run():
        pass
