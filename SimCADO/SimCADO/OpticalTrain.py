###############################################################################
# OpticalTrain 
#
# DESCRIPTION
# The OpticalTrain holds all the information regarding the optical setup as
# well as the individual objects

#
# Classes:
#   
#
# Methods:
#   
#

import sys, os
import warnings

import numpy as np

from astropy.io import fits
import astropy.units as u

import PSFCube as psf
import SpectralCurve as sc
import LightObject as lo
import utils


class OpticalTrain(object):
    
    def __init__(self, cmds):
        self.info = dict([])

        self.cmds = cmds
        
        self.lam_bin_edges   = self.cmds["lam_bin_edges"]
        self.lam_bin_centers = self.cmds["lam_bin_centers"]
        self.lam_res         = self.cmds["SIM_LAM_TC_BIN_WIDTH"]
        self.pix_res         = self.cmds["SIM_INTERNAL_PIX_SCALE"]

    
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
                                                         lam_res=self.lam_res)
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
            fwhm = (1.22*u.rad * self.lam_bin_centers*u.um / \
                                            (m1_diam * u.m)).to(u.arcsec).value 
            if psf_type == "Moffat"
                psf_diff = psf.MoffatPSFCube(self.lam_bin_centers, 
                                             fwhm=fwhm,
                                             pix_res=self.pix_res, 
                                             size=self.psf_size)
            elif psf_type == "Airy"
                psf_diff = psf.AiryPSFCube(self.lam_bin_centers, 
                                           fwhm=fwhm,
                                           pix_res=self.pix_res, 
                                           size=self.psf_size)                                           
                                           
                # Get the Gaussian seeing PSF
                fwhm = (1. - ao_eff/100.) * self.cmds["OBS_SEEING"]
                psf_seeing = psf.GaussianPSFCube(self.lam_bin_centers, 
                                                fwhm=fwhm,
                                                pix_res=self.pix_res)
            
            psf_m1 = psf_diff.convolve(psf_seeing)
        
        scope_psf_master = psf_m1
        
        return scope_psf_master
    
    def make(self, cmds=None):
        """
        To make an optical system, cmds must contain all the keywords from the 
        various .config files. 'cmds' should have been send to __init__, but any
        changes can be passed to make()
        
        Keywords:
        cmds: a dictionary of commands 
        """
    
        self.cmds.update(cmds)
        
        # Make the transmission curve and PSF for the source photons
        self.tc_source  = get_master_tc(self, preset="source")
        self.psf_source = get_master_psf()
        
        # Make the spectral curves for the atmospheric background photons
        self.tc_atmo_bg = get_master_tc(self, preset="atmosphere_bg")
        
        # Get the number of atmospheric background photons in the bandpass
        self.ec_atmo_gb = sc.EmissionCurve(self.cmds["ATMO_EC"])
        
        # Make the transmission curve for the blackbody photons from the mirror
        self.mirror_tc  = get_master_tc(self, preset="mirror_bb")
        
        
        
    
    
    
    
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
    
    
    
        pass
   