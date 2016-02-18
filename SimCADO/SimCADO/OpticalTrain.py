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
import utils


class OpticalTrain(object):
    
    def __init__(self, cmds):
        self.info = dict([])

        self.cmds = cmds
        
        self.lam_bin_edges   = cmds.lam_bin_edges
        self.lam_bin_centers = cmds.lam_bin_centers
        self.lam_res         = cmds.lam_res
        self.pix_res         = cmds.pix_res

        self.lam             = cmds.lam
        self.tc_master       = sc.UnityCurve(lam=self.lam)
        
        self.size            = cmds["SIM_PSF_SIZE"]
        self.psf_master      = psf.DeltaPSFCube(self.lam_bin_centers,
                                                size=self.size,
                                                pix_res=self.pix_res)
        
        self.adc_shifts      = np.zeros((len(self.lam_bin_centers)))
        self.distortion_map  = 1


        
    def read(self, filename):
        pass
    
    def save(self, filename):
        pass
    
    def gen_master_tc(self, tc_keywords=None, preset=None, output=False):
        """
        Combine a list of TransmissionCurves into one, either by specifying the 
        list of command keywords (e.g. ATMO_TC) or by passing a preset keywords
        
        Optional Parameters:
        ===================
        tc_keywords: a list of keywords from the config files. E.g:
                     tc_keywords = ['ATMO_TC', 'SCOPE_M1_TC', 'INST_FILTER_TC']
        preset: a present string for the most common collections of keywords:
                - 'source' includes all the elements seen by source photons
                - 'atmosphere_bg' includes surfaces seen by the atmospheric BG
                - 'mirror_bb' includes surfaces seen by the M1 blackbody photons
        output: [False/True] if True, the master_tc is returned, otherwise, it
                updated the internal parameter self.tc_master
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
            
        tc_master  = sc.UnityCurve()
        for key in tc_list.keys():
            tc_master*= self.tc_list[key]
            
        if output: 
            return tc_master
        else:
            self.tc_master = tc_master

            
    def gen_master_psf(self, psf_type="Airy", output=False):
        """
        Generate a Master PSF for the system. This includes the AO PSF. 
        Notes: Jitter can be applied to detector array as a single PSF, and the 
               ADC shift can be applied to each layer of the PSFCube separately
        
        Parameters
        ==========
        psf_type: 'Moffat', 'Airy'
        output: [False/True] if True, the master_tc is returned, otherwise, it
                updated the internal parameter self.tc_master
        """
       
        ####### PSF CUBES #######
        
        ############################################################
        # !!!!!!!!! USER DEFINED CUBE NOT FULLY TESTED !!!!!!!!!!! #
        # !!!!!!! Analytic still using Airy (x) Gaussian !!!!!!!!! #
        # !!!!!!!!!!!!!!! Implement Moffat PSF !!!!!!!!!!!!!!!!!!! #
        ############################################################
        
        self.psf_size = self.cmds["SIM_PSF_SIZE"]
        
        # Make a PSF for the main mirror. If there is one on file, read it in
        # otherwise generate an Airy+Gaussian (or Moffat, Oliver?)
        
        if self.cmds["SCOPE_USE_PSF_FILE"].lower() != "none" and \
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
                psf_m1 = psf.MoffatPSFCube(self.lam_bin_centers, 
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
        
        psf_master = psf_m1
                   
        if output: 
            return psf_master
        else:
            self.psf_master = psf_master
        
        return scope_psf_master
    
    def gen_adc_shifts(self, output=False):
        """
        
        
        Keywords:
        
        """ 
        para_angle = self.cmds["PARALLACTIC_ANGLE"]
        effectiveness = self.cmds["INST_ADC_EFFICIENCY"] / 100.

        ## get the angle shift for each slice
        angle_shift = [utils.atmospheric_refraction(lam,
                                                    self.cmds["OBS_ZENITH_DIST"],
                                                    self.cmds["ATMO_TEMPERATURE"],
                                                    self.cmds["ATMO_REL_HUMIDITY"],
                                                    self.cmds["ATMO_PRESSURE"],
                                                    self.cmds["SCOPE_LATITUDE"],
                                                    self.cmds["SCOPE_ALTITUDE"])
                       for lam in self.lam_bin_centers]

        ## convert angle shift into number of pixels
        ## pixel shifts are defined with respect to last slice
        pixel_shift = (angle_shift - angle_shift[-1]) / self.pix_res
        if np.max(np.abs(pixel_shift)) > 1000:
            raise ValueError("Pixel shifts too great (>1000), check units")

        ## Rotate by the paralytic angle
        x = -pixel_shift * np.sin(para_angle / 57.29578) * (1. - effectiveness)
        y = -pixel_shift * np.cos(para_angle / 57.29578) * (1. - effectiveness)
        adc_shifts = [(xi, yi) for xi, yi in zip(x, y)]
        
        if output: 
            return adc_shifts
        else:
            self.adc_shifts = adc_shifts

    
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
        self.ec_atmo_gb = sc.EmissionCurve(self.cmds["ATMO_EC"])
        
        # Make the transmission curve for the blackbody photons from the mirror
        self.tc_mirror  = get_master_tc(self, preset="mirror_bb")
        self.ec_mirror  = sc.BlackbodyCurve(lam=self.tc_mirror.lam,
                                            temp=self.cmds["SCOPE_M1_TEMP"])
        
            def __init__(self, lam, temp, **kwargs):
        self.params = { "pix_res" :0.004,
                        "area"    :978,
                        "exptime" :1,
                        "temp"    :273
                        }
        
    
    
    
    
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
        #   - detector  
        
        
        # - the optical path for the atmospheric BG emission [can be exported]
        #   - master transmission curve
        #       - n x mirror
        #       - instrument window
        #       - internal mirrors
        #       - dichroic
        #       - filter
        #       - detector QE

        # - the optical path for the mirror BB emission [can be exported]
        #   - master transmission curve
        #       - n-1 x mirror
        #       - instrument window
        #       - internal mirrors
        #       - dichroic
        #       - filter
        #       - detector QE
class ads        