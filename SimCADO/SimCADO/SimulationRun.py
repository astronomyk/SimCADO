###############################################################################
# SimulationRun
#
# DESCIPTION
#
# File Structure needed to run the simulation
# 
#
#
#
#
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

import numpy as np
from astropy.io import fits
import sys

import PSFCube as psf
import SpectralCurve as sc
import LightObject as lo
import utils


class Simulation(object):
    """
    """
    
    def __init__(self, cmds=None):
        """
        """
        self.cmds = cmds    if cmds is not None     else dict([])
    
    def set_user_commands(self, filename):
        """
        """
        self.cmds = get_user_commands(filename)
        
        
        
    def read_optical_train(self, filename):
        # if we should use a previous optical system, import is
        
    def save_optical_train(self, filename, clobber=False):
        # save the current optical system to a FITS file
        
    def make_optical_train(self):
        # if we should use a previous optical system, import is
        # else generate a new one
        
    def read_source(self, )
        
        
        
    
    
    def run()
        pass
   

def get_user_commands(filename):
    """ 
    1 - get master.config
    2 - read in all the other configs and add them to cmds
    3 - read in the user overrides file
    
    Keywords:
    filename: path to the user.config file
    """
    cmds = utils.read_config("../user_commands/master.config")
    for key in cmds.keys(): 
        cmds.update(utils.read_config(cmds[key]))
    cmds.update(utils.read_config(filename))

    return cmds
   

class OpticalTrain(object):
    
    def __init__(self):
        self.info = dict([])
    
        self.psf_list   = dict([])
        self.tc_list    = dict([])
        self.ec_list    = dict([])
        self.plane_list = dict([])
        self.fpa_list   = dict([])
    
        self.lam_bin_edges   = None
        self.lam_bin_centers = None
        self.lam_res         = None
        self.pix_res         = None

    
    def read(self, filename):
        pass
    
    def save(self, filename):
        pass
    
    def make(self, cmds):
        """
        To make an optical system, cmds must contain all the keywords from the 
        various .config files
        """
    
        # Here we make the optical train. This includes
        # - the optical path for the source photons
        #   - master transmission curve             [can be exported]
        #       - atmosphere
        #       - n x mirror
        #       - instrument window
        #       - internal mirrors
        #       - dichroic
        #       - filter
        #       - detector QE
        #   - list of wave-dep plane effects
        #       - imperfect adc
        #   - master psf cube
        #       - AO psf - analytic or from file    [can be exported]
        #       - jitter psf                        [can be exported]
        #   - list of wave-indep plane effects
        #       - imperfect derotation
        #       - distortion                        [weight map can be exported]
        #       - flat fielding                     [can be exported]
        #
        # - the optical path for the atmospheric BG emission [can be exported]
        #   - master transmission curve
        #       - n x mirror
        #       - instrument window
        #       - internal mirrors
        #       - dichroic
        #       - filter
        #       - detector QE
        #   - list of wave-indep plane effects
        #       - distortion 
        #       - flat fielding
        #       - 
        # - the optical path for the mirror BB emission [can be exported]
        #   - master transmission curve
        #       - n-1 x mirror
        #       - instrument window
        #       - internal mirrors
        #       - dichroic
        #       - filter
        #       - detector QE
        #   - list of wave-indep plane effects
        #       - distortion 
        #       - flat fielding
        #       - 
    
    
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   