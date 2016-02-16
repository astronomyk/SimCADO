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

import PSFCube as psf
import SpectralCurve as sc
import LightObject as lo
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
        # if we should use a previous optical system, import is
        # else generate a new one
        pass
        
    def read_source(self, filename):
        pass
           
    def run():
        pass
        
   
    def set_defaults(self):
        """
        Go through all the input parameters where there could be some ambiguity
        """
        
        # check the output path directory and file name
        if self.cmds["OBS_OUTPUT_DIR"] == "none":
            self.cmds["OBS_OUTPUT_DIR"] = "./"
   
        if self.cmds["OBS_OUTPUT_NAME"] == "none":
            self.cmds["OBS_OUTPUT_NAME"] = "output.fits"
   
        
        # Check if a filter curve file has been give, or a standard broadband name
        if self.cmds["INST_FILTER_TC"] in ["I", "z", "Y", "J", "H", "Ks", "K"]:
            self.cmds["INST_FILTER_TC"] = "../data/TC_filter_" + \
                                            self.cmds["INST_FILTER_TC"] + ".dat"

        # if SIM_USE_FILTER_LAM is true, then use the filter curve to set the
        # wavelength boundaries where the filter is < SIM_FILTER_THRESHOLD
        tc_filt = sc.TransmissionCurve(self.cmds['INST_FILTER_TC'])
        
        if self.cmds["SIM_USE_FILTER_LAM"].lower() == "yes":
            mask = np.where(tc_filt.val > self.cmds["SIM_FILTER_THRESHOLD"])[0]
            lam_min, lam_max = tc_filt.lam[mask[0]], tc_filt.lam[mask[-1]]
            self.cmds["lam_bin_edges"] = np.arange( lam_min, lam_max+1E-7, 
                                                    self.cmds["SIM_LAM_PSF_BIN_WIDTH"])
        else:
            lam_min, lam_max = self.cmds["SIM_LAM_MIN"], self.cmds["SIM_LAM_MAX"]
            self.cmds["lam_bin_edges"] = np.arange( tc_filt.lam[i0], 
                                                    tc_filt.lam[i1]+1E-7, 
                                                    self.cmds["SIM_LAM_TC_BIN_WIDTH"])
        
        self.cmds["lam_bin_centers"] = 0.5 * (self.lam_bin_edges[1:] + \
                                              self.lam_bin_edges[:-1]) 
   
   
   

   
   
   
   
  
   
   
   
   
   
   
   
   