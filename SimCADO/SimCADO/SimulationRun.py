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
        1 - get master.config
        2 - read in all the other configs and add them to cmds
        3 - read in the file with user overrides 
        
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
   
    
   
    def check_defaults(self):
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
   
   
   

class OpticalTrain(object):
    
    def __init__(self, cmds=None):
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

        self.cmds = cmds if cmds is not None else dict([])
    
    def read(self, filename):
        pass
    
    def save(self, filename):
        pass
    
    
    def get_master_tc(self, tc_keywords=None, preset=None):
        """
        Combine a list of TransmissionCurves into one, either by specifying the 
        list of command keywords (e.g. ATMO_TC) or by passing a preset keywords
        
        Keywords:
        tc_keywords: a list of keywords from the config files. E.g:
                     tc_keywords = ['ATMO_TC', 'SCOPE_M1_TC', 'INST_FILTER_TC']
        preset: a present string for the most common collections of keywords:
                - 'source' includes all the elements seen by source photons
                - 'atmosphere_bg' includes surfaces seen by the atmospheric BG
                - 'mirror_bb' includes surfaces seen by the M1 blackbody photons
        """
        
        if tc_keywords is None:
            if preset == "source":
                tc_keywords = ['SCOPE_M1_TC']*self.cmds['SCOPE_NUM_MIRRORS'] + \
                              ['ATMO_TC', 'INST_ADC_TC', 'INST_DICHROIC_TC', 
                               'INST_ENTR_WINDOW_TC', 'INST_FILTER_TC', 'FPA_QE']
            if preset == "atmosphere_bg":
                tc_keywords = ['SCOPE_M1_TC']*self.cmds['SCOPE_NUM_MIRRORS'] + \
                              ['INST_ADC_TC', 'INST_DICHROIC_TC', 
                               'INST_ENTR_WINDOW_TC', 'INST_FILTER_TC', 'FPA_QE']
            if preset == "mirror_bb":
                tc_keywords = ['SCOPE_M1_TC']*(self.cmds['SCOPE_NUM_MIRRORS']-1) + \
                              ['INST_ADC_TC', 'INST_DICHROIC_TC', 
                               'INST_ENTR_WINDOW_TC', 'INST_FILTER_TC', 'FPA_QE']
        
        for key in tc_keywords:
            if key not in self.cmds.keys():
                raise ValueError(key+" is not in your list of commands")
            
            if self.cmds[key].lower() != 'none':     
                self.tc_list[key] = sc.TransmissionCurve(filename=self.cmds[key],
                                            lam_res=cmds["SIM_LAM_TC_BIN_WIDTH"])
            else: 
                self.tc_list[key] = sc.UnityCurve()
            
        tc_master_scope  = sc.UnityCurve()
        for key in tc_list.keys():
            tc_master_scope *= self.tc_list[key]
            
        return tc_master_scope

    def get_master_psf(self, psf_type="Airy"):
        """
        Generate a Master PSF for the system. This includes the AO PSF. 
        Notes: Jitter can be applied to detector array as a single PSF, and the 
               ADC shift can be applied to each layer of the PSFCube seperately
        
        Parameters
        ==========
        psf_type: 'Moffat', 'Airy'
        """
       
        ####### PSF CUBES #######
        
        ############################################################
        # !!!!!!!!! USER DEFINED CUBE NOT FULLY TESTED !!!!!!!!!!! #
        # !!!!!!! Analytic still using Airy (x) Gaussian !!!!!!!!! #
        # !!!!!!!!!!!!!!! Implement Moffat PSF !!!!!!!!!!!!!!!!!!! #
        ############################################################
        
        self.pix_res = self.cmds["SIM_INTERNAL_PIX_SCALE"]
        self.psf_size = self.cmds["SIM_PSF_SIZE"]
        self.area = np.pi * (self.cmds["SCOPE_M1_DIAMETER_OUT"]**2 - \
                             self.cmds["SCOPE_M1_DIAMETER_IN"]**2)
        
        # Make a PSF for the main mirror. If there is one on file, read it in
        # otherwise generate an Airy+Gaussian (or Moffat, Oliver?)
        
        if self.cmds["SCOPE_USE_PSF_FILE"] != "none" and \
                                os.path.exists(self.cmds["SCOPE_USE_PSF_FILE"]):
            psf_m1 = psf.UserPSFCube(self.cmds["SCOPE_USE_PSF_FILE"])
            if psf_m1[0].pix_res != self.pix_res:
                psf_m1 = psf_m1.resample(self.pix_res)    
        else:
            m1_diam = self.cmds["SCOPE_M1_DIAMETER_OUT"]
            ao_eff  = self.cmds["SCOPE_AO_EFFECTIVENESS"]
            
            # Get a Diffraction limited PSF
            fwhm = (1.22*u.rad * self.cmds["lam_bin_centers"]*u.um / \
                                            (m1_diam * u.m)).to(u.arcsec).value 
            if psf_type == "Moffat"
                psf_diff = psf.MoffatPSFCube(self.cmds["lam_bin_centers"], 
                                             fwhm=fwhm,
                                             pix_res=self.pix_res, 
                                             size=self.psf_size)
            elif psf_type == "Airy"
                psf_diff = psf.AiryPSFCube(self.cmds["lam_bin_centers"], 
                                           fwhm=fwhm,
                                           pix_res=self.pix_res, 
                                           size=self.psf_size)                                           
                                           
                # Get the Gaussian seeing PSF
                fwhm = (1. - ao_eff/100.) * self.cmds["OBS_SEEING"]
                psf_seeing = psf.GaussianPSFCube(self.cmds["lam_bin_centers"], 
                                                fwhm=fwhm,
                                                pix_res=self.pix_res)
            
            psf_m1 = psf_diff.convolve(psf_seeing)
        
        scope_psf_master = psf_m1
        
        return scope_psf_master
    
    def make(self, cmds=None):
        """
        To make an optical system, cmds must contain all the keywords from the 
        various .config files
        """
    
        self.cmds.update(cmds)
    
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


            
        
        
        
        
        
        ####################################
        # Idea - break these up into make_tc(), make_psf() and generic make()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
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
    
    
        SCOPE_M1_TEMP
    
    
    
        pass
   
   
   
   
   
  
   
   
   
   
   
   
   
   