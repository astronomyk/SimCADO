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
try:
    import SimCADO.Detector as fpa
    import SimCADO.PSFCube as psf
    import SimCADO.SpectralCurve as sc
    import SimCADO.PlaneEffect as pe
    import SimCADO.utils as utils
except:
    import Detector as fpa
    import PSFCube as psf
    import SpectralCurve as sc
    import PlaneEffect as pe
    import utils

class OpticalTrain(object):

    def __init__(self, cmds):
        self.info = dict([])

        self.cmds = cmds

        fname = self.cmds["OPTICAL_TRAIN_IN_PATH"]
        if fname is not None:
            if not os.path.exists(fname):
                raise ValueError(fname+" doesn't exist")
            
            self.read(fname)
        else:
            
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

            self.adc_shifts      = np.zeros((len(self.lam_bin_centers),2))
            #self.distortion_map

            self.make()
        
        
    def make(self, cmds=None):
        """
        To make an optical system, cmds must contain all the keywords from the
        various .config files. 'cmds' should have been send to __init__, but any
        changes can be passed to make()

        Parameters
        ==========
        - cmds : UserCommands
            a dictionary of commands
        """

        self.cmds.update(cmds)

        # Make the transmission curve for the blackbody photons from the mirror
        self.tc_mirror  = self._gen_master_tc(preset="mirror_bb")
        self.ec_mirror  = sc.BlackbodyCurve(lam=self.tc_mirror.lam,
                                            temp=self.cmds["SCOPE_M1_TEMP"])

        # Make the spectral curves for the atmospheric background photons
        self.tc_atmo_bg = self._gen_master_tc(preset="atmosphere_bg")
        self.ec_atmo_gb = sc.EmissionCurve(filename = self.cmds["ATMO_EC"])

        # Make the transmission curve and PSF for the source photons
        self.tc_source  = self._gen_master_tc(preset="source")
        self.psf_source = self._gen_master_psf()

        # Make a detector Plane
        self.detector = self._gen_detector()
        
        # Get the ADC shifts, telescope shake and field rotation angle
        self.adc_shifts = self._gen_adc_shifts()
        self.jitter_psf = self._gen_telescope_shake(self):
        # self.field_rot = self._gen_field_rotation_angle()
        
        

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

    def read(self, filename):
        pass

    def save(self, filename):
        pass
        
        

    def _gen_master_tc(self, tc_keywords=None, preset=None):
        """
        Combine a list of TransmissionCurves into one, either by specifying the
        list of command keywords (e.g. ATMO_TC) or by passing a preset keywords

        Optional Parameters:
        ===================
        tc_keywords: a list of keywords from the .config files. E.g:
                     tc_keywords = ['ATMO_TC', 'SCOPE_M1_TC', 'INST_FILTER_TC']
        preset: a present string for the most common collections of keywords:
                - 'source' includes all the elements seen by source photons
                - 'atmosphere_bg' includes surfaces seen by the atmospheric BG
                - 'mirror_bb' includes surfaces seen by the M1 blackbody photons
        """

        if tc_keywords is None:
            if preset is not None:
                base = ['SCOPE_M1_TC'] * (int(self.cmds['SCOPE_NUM_MIRRORS']) - 1) + \
                       ['INST_ADC_TC', 'INST_DICHROIC_TC', 'INST_ENTR_WINDOW_TC', 
                        'INST_FILTER_TC', 'FPA_QE']
                if preset == "source":
                    tc_keywords = ['ATMO_TC'] + ['SCOPE_M1_TC'] + base
                if preset == "atmosphere_bg":
                    tc_keywords = ['SCOPE_M1_TC'] + base
                if preset == "mirror_bb":
                    tc_keywords = base
            else:
                warnings.warn("""
                No presets or keywords passed to gen_master_tc(). 
                Setting self.tc_master = sc.UnityCurve()""")
                self.tc_master = sc.UnityCurve()
                return
                
        tc_dict = dict([])
                
        for key in tc_keywords:
            if key not in self.cmds.keys():
                raise ValueError(key + " is not in your list of commands")

            if self.cmds[key].lower() != 'none':
                tc_dict[key] = sc.TransmissionCurve(filename=self.cmds[key],
                                                    lam_res=self.lam_res)
            else:
                tc_dict[key] = sc.UnityCurve()

        tc_master = sc.UnityCurve( lam=self.lam, lam_res=self.lam_res, 
                                   min_step=self.cmds["SIM_SPEC_MIN_STEP"])
        for key in tc_keywords:
            tc_master *= tc_dict[key]

        return tc_master


    def _gen_master_psf(self, psf_type="Airy"):
        """
        Generate a Master PSF for the system. This includes the AO PSF.
        Notes: Jitter can be applied to detector array as a single PSF, and the
               ADC shift can be applied to each layer of the PSFCube separately

        Parameters
        ==========
        psf_type : str
            'Moffat', 'Airy'
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

        if self.cmds["SCOPE_PSF_FILE"].lower() != "none" and \
                                os.path.exists(self.cmds["SCOPE_PSF_FILE"]):
            psf_m1 = psf.UserPSFCube(self.cmds["SCOPE_PSF_FILE"])
            if psf_m1[0].pix_res != self.pix_res:
                psf_m1 = psf_m1.resample(self.pix_res)
        else:
            m1_diam = self.cmds["SCOPE_M1_DIAMETER_OUT"]
            ao_eff  = self.cmds["SCOPE_AO_EFFECTIVENESS"]

            # Get a Diffraction limited PSF
            fwhm = (1.22*u.rad * self.lam_bin_centers * u.um / \
                                            (m1_diam * u.m)).to(u.arcsec).value
            if psf_type == "Moffat":
                psf_m1 = psf.MoffatPSFCube(self.lam_bin_centers,
                                           fwhm=fwhm,
                                           pix_res=self.pix_res,
                                           size=self.psf_size)
            elif psf_type == "Airy":
                psf_diff = psf.AiryPSFCube(self.lam_bin_centers,
                                           fwhm=fwhm,
                                           pix_res=self.pix_res,
                                           size=self.psf_size)

                # Get the Gaussian seeing PSF
                if ao_eff < 100.:
                    fwhm = (1. - ao_eff/100.) * self.cmds["OBS_SEEING"]
                    psf_seeing = psf.GaussianPSFCube(self.lam_bin_centers,
                                                    fwhm=fwhm,
                                                    pix_res=self.pix_res)

                    psf_m1 = psf_diff.convolve(psf_seeing)
                else:
                    psf_m1 = psf_diff
                    
        psf_master = psf_m1
        return psf_master


    def _gen_detector(self):
        """
        Keywords:
        """
        fname = self.cmds["FPA_NOISE_PATH"]
        if fname is not None:
            if not os.path.exists(fname):
                raise ValueError(fname+" doesn't exist")
        
        return fpa.Detector(fname, **self.cmds)
        

    def _gen_adc_shifts(self):
        """
        Keywords:
        """
        adc_shifts = pe.adc_shift(self.cmds)       
        return adc_shifts


    def _gen_field_rotation_angle(self):
        pass
        
    def _gen_telescope_shake(self):
        """
        Keywords:
        """
        pix_res =   self.cmds["SIM_DETECTOR_PIX_SCALE"] / \
                    self.cmds["SIM_OVERSAMPLING"]
        jitter_psf = psf.GaussianPSF(   fwhm=self.cmds["SCOPE_JITTER_FWHM"], 
                                        pix_res=pix_res)
        return jitter_psf
    

class ads:
    pass
