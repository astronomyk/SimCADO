"""optics.py"""
###############################################################################
# OpticalTrain
#
# DESCRIPTION
# The OpticalTrain holds all the information regarding the optical path as
# well as the individual objects
#
# TODO List
# =========
# - Make the Detector independent of the OpticalTrain
# - Implement saving and reloading of OpticalTrain objects
#

import os

import glob
import warnings, logging
from datetime import datetime as dt
from copy import deepcopy

import numpy as np

from astropy.io import fits
from astropy.io import ascii as ioascii    # 'ascii' redefines built-in
import astropy.units as u

from . import psf as psf
from . import spectral as sc
from . import spatial as pe
from .commands import UserCommands
from .utils import __pkg_dir__

__all__ = ["OpticalTrain", "get_filter_curve", "get_filter_set"]

class OpticalTrain(object):
    """
    The OpticalTrain object reads in or generates the information necessary to
    model the optical path for all (3) sources of photons: the astronomical
    source, the atmosphere and the primary mirror.

    Parameters
    ----------
    cmds : UserCommands, optional
        Holds the commands needed to generate a model of the optical train


    Optional Parameters
    -------------------
    Any keyword-value pair contained in the default configuration file


    See Also
    --------
    .commands.dump_defaults(), .commands.UserCommands


    General Attributes
    ------------------
    - cmds : commands, optional
        a dictionary of commands for running the simulation

    Spatial attributes (for PSFs)
    -----------------------------
    lam_bin_edges : 1D array
        [um] wavelengths of the edges of the bins used for each PSF layer
    lam_bin_centers
        [um] wavelengths of the centre of the bins used for each PSF layer
    pix_res : float
        [arcsec] infernal oversampled pixel resolution (NOT detector plate scale)
    psf :
    jitter_psf :
    adc_shifts
        [pixel]
    field_rot :
        [degrees]

    Spectral attributes (for EmissionCurves)
    ----------------------------------------
    lam : 1D array
        [um] Vector of wavelength bins for spectrals
    lam_res : float
        [um] resolution between
    psf_size : int
        [pixels] The width of a PSF
    tc_ao : TransmissionCurve
        [0..1]
    tc_mirror : TransmissionCurve
        [0..1]
    tc_atmo : TransmissionCurve
        [0..1]
    tc_source : TransmissionCurve
        [0..1]
    ec_ao : EmissionCurve
        [ph/s/voxel]
    ec_mirror : EmissionCurve
        [ph/s/voxel]
    ec_atmo : EmissionCurve
        [ph/s/voxel]
    ph_ao : EmissionCurve
        [ph/s/voxel]
    ph_mirror : EmissionCurve
        [ph/s/voxel]
    ph_atmo : EmissionCurve
        [ph/s/voxel]
    n_ph_ao : float
        [ph/s]
    n_ph_mirror : float
        [ph/s]
    n_ph_atmo : float
        [ph/s]


    """

    def __init__(self, cmds, **kwargs):
        self.info = dict([])

        if cmds is None:
            cmds = UserCommands()
        self.cmds = deepcopy(cmds)
        self.cmds.update(kwargs)

        self.tc_master = None   # set in separate method
        self.psf_size = None   # set in separate method

        fname = self.cmds["SIM_OPT_TRAIN_IN_PATH"]
        if fname is not None:
            if not os.path.exists(fname):
                raise ValueError(fname+" doesn't exist")

            self.read(fname)
        else:
            self.lam_bin_edges = cmds.lam_bin_edges
            self.lam_bin_centers = cmds.lam_bin_centers
            self.pix_res = cmds.pix_res

            self.lam = cmds.lam
            self.lam_res = cmds.lam_res

            self._make()


    def _make(self, cmds=None):
        """
        To make an optical system, cmds must contain all the keywords from the
        various .config files. 'cmds' should have been send to __init__, but any
        changes can be passed to make()

        Parameters
        ----------
        - cmds : commands
            a dictionary of commands
        """

        # Here we make the optical train. This includes
        # - the optical path for the source photons
        #   - master transmission curve             [can be exported]
        #       - atmosphere, n x mirror, instrument window, internal mirrors
        #       - dichroic, filter, detector QE
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
        #       - noise frame,

        if self.cmds.verbose:
            print("Generating an optical train")
        logging.debug("[OpticalTrain] Generating an optical train")

        if cmds is not None:
            self.cmds.update(cmds)
        self._load_all_tc()
        self._gen_all_tc()
        self.psf = self._gen_master_psf()

        # Get the ADC shifts, telescope shake and field rotation angle
        self.adc_shifts = self._gen_adc_shifts()
        self.jitter_psf = self._gen_telescope_shake()
        #self.field_rot = self._gen_field_rotation_angle()


    def replace_psf(self, psf, lam_bin_centres):
        """
        Change the PSF of the optical train
        """
        pass


    def update_filter(self, trans=None, lam=None, filter_name=None):
        """
        Update the filter curve without recreating the full OpticalTrain object

        Parameters
        ----------
        trans : TransmissionCurve, np.array, list, optional
            [0 .. 1] the transmission coefficients. Either a TransmissionCurve
            object can be passed (in which case omit ``lam``) or an array/list can
            be passed (in which case specify ``lam``)
        lam : np.array, list, optional
            [um] an array for the spectral bin centres, if ``trans`` is not a
            TransmissionCurve object
        filter_name : str, optional
            The name of a filter curve contained in the package_dir. User
            get_filter_set() to find which filter curves are installed.

        See also
        --------
        :class:`simcado.spectral.TransmissionCurve`
        :func:`simcado.optics.get_filter_set`

        """
        if filter_name == lam == trans == None:
            raise ValueError("At least one parameter must be specified")

        if filter_name is not None:
            filt = get_filter_curve(filter_name)
        elif trans is not None:
            if isinstance(trans, (sc.TransmissionCurve,
                                  sc.EmissionCurve,
                                  sc.UnityCurve,
                                  sc.BlackbodyCurve)):
                filt = trans
            elif isinstance(trans, (np.ndarray, list, tuple)) and \
                 isinstance(lam  , (np.ndarray, list, tuple)):
                filt = sc.TransmissionCurve(lam=lam, val=trans,
                                            lam_res=self.lam_res)

        self.cmds["INST_FILTER_TC"] = filt
        self._gen_all_tc()


    def read(self, filename):
        pass

    def save(self, filename):
        pass


    def apply_tracking(self, arr):
        return pe.tracking(arr, self.cmds)

    def apply_derotator(self, arr):
        return pe.derotator(arr, self.cmds)

    def apply_wind_jitter(self, arr):
        return pe.wind_jitter(arr, self.cmds)

    def _load_all_tc(self, tc_list=["ATMO_TC", "SCOPE_M1_TC", "INST_MIRROR_AO_TC",
                                    "INST_ENTR_WINDOW_TC", "INST_DICHROIC_TC",
                                    "INST_MIRROR_TC", "INST_ADC_TC",
                                    "INST_PUPIL_TC", "INST_FILTER_TC", "FPA_QE"]):
        """
        Pre-loads all the transmission curves
        """
        for tc in tc_list:
            if isinstance(self.cmds[tc], str):
                airmass = self.cmds["ATMO_AIRMASS"] if tc == "ATMO_TC" else None
                self.cmds[tc] = sc.TransmissionCurve(filename=self.cmds[tc],
                                                     airmass=airmass)
            elif self.cmds[tc] is None:
                self.cmds[tc] = sc.UnityCurve()

        # see Rics email from 22.11.2016
        wfe = self.cmds["INST_TOTAL_WFE"]
        lam = np.arange(0.3,3.0)
        val = np.exp( -(2 * np.pi * (wfe*u.nm) / (lam*u.um))**2 )
        self.cmds.cmds["INST_SURFACE_FACTOR"] = sc.TransmissionCurve(lam=lam,
                                                                     val=val)


    def _gen_all_tc(self):

        ############## AO INSTRUMENT PHOTONS #########################
        if self.cmds.verbose:
            print("Generating AO module mirror emission photons")
        logging.debug("[_gen_all_tc] Generating AO module mirror emission photons")

        # get the total area of mirrors in the telescope
        # !!!!!! Bad practice, this is E-ELT specific hard-coding !!!!!!
        
        mirr_list = self.cmds.mirrors_ao
        ao_area = np.pi / 4 * np.sum(mirr_list["Outer"]**2 - \
                                     mirr_list["Inner"]**2)

        # Make the transmission curve for the blackbody photons from the mirror
        self.tc_ao = self._gen_master_tc(preset="ao")
        self.ec_ao = sc.BlackbodyCurve(lam     = self.tc_ao.lam,
                                       temp    = self.cmds["INST_AO_TEMPERATURE"],
                                       pix_res = self.cmds.pix_res,
                                       area    = ao_area)

        if self.cmds["INST_USE_AO_MIRROR_BG"].lower() == "yes" and \
           self.cmds["SCOPE_PSF_FILE"].lower() != "scao":

            self.ph_ao = self.ec_ao * self.tc_ao
            self.n_ph_ao = self.ph_ao.photons_in_range(self.lam_bin_edges[0],
                                                       self.lam_bin_edges[-1])
        else:
            self.ec_ao = None
            self.ph_ao = None
            self.n_ph_ao = 0.



        ############## TELESCOPE PHOTONS #########################
        if self.cmds.verbose:
            print("Generating telescope mirror emission photons")
        logging.debug("[_gen_all_tc] Generating telescope mirror emission photons")

        # get the total area of mirrors in the telescope
        # !!!!!! Bad practice, this is E-ELT specific hard-coding !!!!!!
        mirr_list = self.cmds.mirrors_telescope
        scope_area = np.pi / 4 * np.sum(mirr_list["Outer"]**2 - \
                                        mirr_list["Inner"]**2)

        # Make the transmission curve for the blackbody photons from the mirror
        self.tc_mirror = self._gen_master_tc(preset="mirror")
        self.ec_mirror = sc.BlackbodyCurve(lam     = self.tc_mirror.lam,
                                           temp    = self.cmds["SCOPE_TEMP"],
                                           pix_res = self.cmds.pix_res,
                                           area    = scope_area)

        if self.cmds["SCOPE_USE_MIRROR_BG"].lower() == "yes":
            self.ph_mirror = self.ec_mirror * self.tc_mirror
            self.n_ph_mirror = self.ph_mirror.photons_in_range(self.lam_bin_edges[0],
                                                               self.lam_bin_edges[-1])
        else:
            self.ec_mirror = None
            self.ph_mirror = None
            self.n_ph_mirror = 0.


        ############## ATMOSPHERIC PHOTONS #########################
        if self.cmds.verbose:
            print("Generating atmospheric emission photons")
        logging.debug("[_gen_all_tc] Generating atmospheric emission photons")

        # Make the spectral curves for the atmospheric background photons
        self.tc_atmo = self._gen_master_tc(preset="atmosphere")

        if self.cmds["ATMO_USE_ATMO_BG"].lower() == "yes":

            if self.cmds["ATMO_EC"] is not None:
                self.ec_atmo = sc.EmissionCurve(filename=self.cmds["ATMO_EC"],
                                                pix_res=self.cmds.pix_res,
                                                area=self.cmds.area,
                                                airmass=self.cmds["ATMO_AIRMASS"])
                self.ph_atmo = self.tc_atmo * self.ec_atmo
                self.n_ph_atmo = self.ph_atmo.photons_in_range(self.lam_bin_edges[0],
                                                               self.lam_bin_edges[-1])
            else:
                if len(self.cmds["INST_FILTER_TC"]) > 2:
                    filt = self.cmds["INST_FILTER_TC"][-5:-3].replace(".", "")
                else:
                    filt = self.cmds["INST_FILTER_TC"]

                if filt not in "BVRIzYJHKKs":
                    raise ValueError("""Only broadband filters (BVRIzYJHKKs)
                      can be used with keyword ATMO_BG_MAGNITUDE. Please provide
                      a filename for ATMO_EC, or write ATMO_EC  default""")

                if self.cmds["ATMO_BG_MAGNITUDE"] == "default":
                    path = os.path.join(__pkg_dir__, "data", "EC_sky_magnitude.dat")
                    self.cmds["ATMO_BG_MAGNITUDE"] = ioascii.read(path)[filt][0]


                from simcado import source
                ph_zero = source.zero_magnitude_photon_flux(filt)

                self.ec_atmo = None
                self.ph_atmo = None

                factor = 10**(-0.4*self.cmds["ATMO_BG_MAGNITUDE"]) * \
                                        self.cmds.pix_res**2 * self.cmds.area

                self.n_ph_atmo = ph_zero * factor

        else:
            self.ec_atmo = None
            self.ph_atmo = None
            self.n_ph_atmo = 0.

        self.n_ph_bg = self.n_ph_atmo + self.n_ph_mirror +  self.n_ph_ao

        ############## SOURCE PHOTONS #########################
        if self.cmds.verbose:
            print("Generating optical path for source photons")
        logging.debug("[_gen_all_tc] GGenerating optical path for source photons")

        # Make the transmission curve and PSF for the source photons
        self.tc_source = self._gen_master_tc(preset="source")


    def _gen_master_tc(self, tc_keywords=None, preset=None):
        """
        Combine a list of TransmissionCurves into one, either by specifying the
        list of command keywords (e.g. ATMO_TC) or by passing a preset keywords

        Optional Parameters
        ===================
        tc_keywords: a list of keywords from the .config files. E.g:
                     tc_keywords = ['ATMO_TC', 'SCOPE_M1_TC', 'INST_FILTER_TC']
        preset: a present string for the most common collections of keywords:
                - 'source' includes all the elements seen by source photons
                - 'atmosphere' includes surfaces seen by the atmospheric BG
                - 'mirror' includes surfaces seen by the M1 blackbody photons
                - 'ao' inludes surfaces seen by the AO module blackbody photons
        """

        if tc_keywords is None:
            if preset is not None:
                base = ['INST_MIRROR_AO_TC'] * (int(self.cmds['INST_NUM_AO_MIRRORS']) - 1) + \
                       ['INST_ENTR_WINDOW_TC'] * int(self.cmds['INST_ENTR_NUM_SURFACES']) + \
                       ['INST_DICHROIC_TC']  * int(self.cmds['INST_DICHROIC_NUM_SURFACES']) + \
                       ['INST_MIRROR_TC']    * (int(self.cmds['INST_NUM_MIRRORS'])) + \
                       ['INST_ADC_TC']       * int(self.cmds['INST_ADC_NUM_SURFACES']) + \
                       ['INST_PUPIL_TC']     * int(self.cmds['INST_PUPIL_NUM_SURFACES']) + \
                       ['INST_FILTER_TC'] + \
                       ['INST_SURFACE_FACTOR'] + \
                       ['FPA_QE']

                ao =   ['INST_MIRROR_AO_TC'] + ['SCOPE_M1_TC'] * (int(self.cmds['SCOPE_NUM_MIRRORS']) - 1)


                if preset == "ao":
                    tc_keywords = base
                elif preset == "mirror":
                    tc_keywords = base + ao
                elif preset == "atmosphere":
                    tc_keywords = base + ao + ['SCOPE_M1_TC']
                elif preset == "source":
                    tc_keywords = base + ao + ['SCOPE_M1_TC', 'ATMO_TC']
                else:
                    raise ValueError("Unknown preset parameter " + preset)

            else:
                warnings.warn("""
                No presets or keywords passed to gen_master_tc().
                Setting self.tc_master = sc.UnityCurve()""", UserWarning)
                self.tc_master = sc.UnityCurve()
                return

        tc_dict = dict([])

        for key in tc_keywords:
            if key not in self.cmds.keys():
                raise ValueError(key + " is not in your list of commands")

            if self.cmds[key] is not None:
                if isinstance(self.cmds[key], (sc.TransmissionCurve,
                                               sc.EmissionCurve,
                                               sc.UnityCurve,
                                               sc.BlackbodyCurve)):
                    tc_dict[key] = self.cmds[key]
                else:
                    airmass = self.cmds["ATMO_AIRMASS"] if key == "ATMO_TC" else None
                    tc_dict[key] = sc.TransmissionCurve(filename=self.cmds[key],
                                                        lam_res=self.lam_res,
                                                        airmass=airmass)
            else:
                tc_dict[key] = sc.UnityCurve()

        tc_master = sc.UnityCurve(lam=self.lam, lam_res=self.lam_res,
                                  min_step=self.cmds["SIM_SPEC_MIN_STEP"])
        for key in tc_keywords:
            tc_master *= tc_dict[key]

        return tc_master


    def _gen_master_psf(self):
        """
        Import or make a aaster PSF for the system. 
        
        Notes
        -----
        Jitter can be applied to detector array as a single PSF, and the
               ADC shift can be applied to each layer of the psf separately

        """

        ############################################################
        # !!!!!!!!! USER DEFINED CUBE NOT FULLY TESTED !!!!!!!!!!! #
        # !!!!!!! Analytic still using Airy (x) Gaussian !!!!!!!!! #
        # !!!!!!!!!!!!!!! Implement Moffat PSF !!!!!!!!!!!!!!!!!!! #
        ############################################################

        self.psf_size = self.cmds["SIM_PSF_SIZE"]

        # Make a PSF for the main mirror. If there is one on file, read it in
        # otherwise generate an Airy+Gaussian (or Moffat, Oliver?)
        
        if self.cmds["SCOPE_PSF_FILE"] is None:
            warnings.warn("""
            No PSF given. SCOPE_PSF_FILE = None. 
            Returning an Delta function for SCOPE_PSF_FILE""")

            psf_m1 = psf.DeltaPSFCube(self.lam_bin_centers,
                                      pix_res=self.pix_res,
                                      size=9)
            logging.debug("No PSF Given: making Delta PSF")
            
        elif isinstance(self.cmds["SCOPE_PSF_FILE"], psf.PSFCube):
            psf_m1 = self.cmds["SCOPE_PSF_FILE"]
            logging.debug("Using PSF: " + self.cmds["SCOPE_PSF_FILE"])
        
        elif isinstance(self.cmds["SCOPE_PSF_FILE"], str):
            if self.cmds.verbose:
                print("Using PSF:", self.cmds["SCOPE_PSF_FILE"])
            
            if os.path.exists(self.cmds["SCOPE_PSF_FILE"]):
                logging.debug("Using PSF: " + self.cmds["SCOPE_PSF_FILE"])
            
                psf_m1 = psf.UserPSFCube(self.cmds["SCOPE_PSF_FILE"],
                                         self.lam_bin_centers)

                if psf_m1[0].pix_res != self.pix_res:
                    psf_m1.resample(self.pix_res)
            else:
                warnings.warn("""
                Couldn't resolve SCOPE_PSF_FILE. 
                Returning an Delta function for SCOPE_PSF_FILE""")
        
                psf_m1 = psf.DeltaPSFCube(self.lam_bin_centers,
                                          pix_res=self.pix_res,
                                          size=9)
                logging.debug("Couldn't resolve given PSF: making Delta PSF")
            
        return psf_m1

        
    def _gen_adc_shifts(self):
        """
        Keywords:
        """
        adc_shifts = pe.adc_shift(self.cmds)
        return adc_shifts


    def _gen_field_rotation_angle(self):
        return 0


    def _gen_telescope_shake(self):
        """
        Keywords:
        """
        jitter_psf = psf.GaussianPSF(fwhm=self.cmds["SCOPE_JITTER_FWHM"],
                                     pix_res=self.cmds.pix_res)
        return jitter_psf


## note: 'filter' redefines a built-in and should not be used
def get_filter_curve(filter_name):
    """
    Return a Vis/NIR broadband filter TransmissionCurve object

    Parameters
    ----------
    filter_name : str

    Notes
    -----
    Acceptable filters can be found be calling get_filter_set()
    """

    if filter_name not in get_filter_set(path=None):
        raise ValueError("filter not recognised: "+filter)
    fname = os.path.join(__pkg_dir__, "data", "TC_filter_"+filter_name+".dat")
    return sc.TransmissionCurve(filename=fname)


def get_filter_set(path=None):
    """
    Return a list of the filters installed in the package directory
    """
    if path is None:
        path = os.path.join(__pkg_dir__, "data")
    lst = [i.replace(".dat", "").split("TC_filter_")[-1] \
                    for i in glob.glob(os.path.join(path, "TC_filter*.dat"))]
    return lst
