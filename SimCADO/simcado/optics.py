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
import inspect
__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))

import glob
import warnings
from datetime import datetime as dt

import numpy as np

from astropy.io import fits
from astropy.io import ascii as ioascii    # 'ascii' redefines built-in
import astropy.units as u

try:
    #import simcado.detector as fpa
    import simcado.psf as psf
    import simcado.spectral as sc
    import simcado.spatial as pe
    #import simcado.utils as utils
except ImportError:
    #import detector as fpa
    import psf as psf
    import spectral as sc
    import spatial as pe
    #import utils

__all__ = ["OpticalTrain"]

class OpticalTrain(object):
    """
    The OpticalTrain object reads in or generates the information necessary to
    model the optical path for all (3) sources of photons: the astronomical
    source, the atmosphere and the primary mirror.

    Attributes
    ==========
    General attributes
    ------------------
    - cmds : commands
        a dictionary of commands for running the simulation

    Spatial attributes (for PSFs)
    -----------------------------
    - lam_bin_edges : 1D array
        [um] wavelengths of the edges of the bins used for each PSF layer
    - lam_bin_centers
        [um] wavelengths of the centre of the bins used for each PSF layer
    - pix_res : float
        [arcsec] infernal oversampled pixel resolution (NOT detector plate scale)
    - psf :
    - jitter_psf :
    - adc_shifts
        [pixel]
    - field_rot :
        [degrees]

    Spectral attributes (for EmissionCurves)
    ----------------------------------------
    - lam : 1D array
        [um] Vector of wavelength bins for spectrals
    - lam_res : float
        [um] resolution between
    - psf_size : int
        [pixels] The width of a PSF
    - tc_ao : TransmissionCurve
        [0..1]
    - tc_mirror : TransmissionCurve
        [0..1]
    - tc_atmo : TransmissionCurve
        [0..1]
    - tc_source : TransmissionCurve
        [0..1]
    - ec_ao : EmissionCurve
        [ph/s/voxel]
    - ec_mirror : EmissionCurve
        [ph/s/voxel]
    - ec_atmo : EmissionCurve
        [ph/s/voxel]
    - ph_ao : EmissionCurve
        [ph/s/voxel]
    - ph_mirror : EmissionCurve
        [ph/s/voxel]
    - ph_atmo : EmissionCurve
        [ph/s/voxel]
    - n_ph_ao : float
        [ph/s]
    - n_ph_mirror : float
        [ph/s]
    - n_ph_atmo : float
        [ph/s]


    """
    def __init__(self, cmds):
        self.info = dict([])

        self.cmds = cmds

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
        if cmds is not None:
            self.cmds.update(cmds)


        ############## AO INSTRUMENT PHOTONS #########################    
        if self.cmds.verbose:
            print("Generating AO module mirror emission photons")

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
        # Make the spectral curves for the atmospheric background photons
        self.tc_atmo = self._gen_master_tc(preset="atmosphere")

        if self.cmds["ATMO_USE_ATMO_BG"].lower() == "yes":

            if self.cmds["ATMO_EC"] is not None:
                self.ec_atmo = sc.EmissionCurve(filename=self.cmds["ATMO_EC"],
                                                pix_res=self.cmds.pix_res,
                                                area=self.cmds.area)
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
        # Make the transmission curve and PSF for the source photons
        self.tc_source = self._gen_master_tc(preset="source")
        self.psf = self._gen_master_psf()

        # Detector has been outsourced - this is irrelevant now
          ############## detector #########################
        #if self.cmds.verbose:
        #    print("Generating the detector array")
        # Make a detector Plane
        #self.detector = fpa.Detector(self.cmds)
        #self.chips = self.detector.chips

        # Get the ADC shifts, telescope shake and field rotation angle
        self.adc_shifts = self._gen_adc_shifts()
        self.jitter_psf = self._gen_telescope_shake()
        #self.field_rot = self._gen_field_rotation_angle()


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
                if preset == "mirror":
                    tc_keywords = base + ao
                if preset == "atmosphere":
                    tc_keywords = base + ao + ['SCOPE_M1_TC']
                if preset == "source":
                    tc_keywords = base + ao + ['SCOPE_M1_TC', 'ATMO_TC']
                
            else:
                raise ValueError("""
                No presets or keywords passed to gen_master_tc().
                Setting self.tc_master = sc.UnityCurve()""")
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
                    tc_dict[key] = sc.TransmissionCurve(filename=self.cmds[key],
                                                    lam_res=self.lam_res)
            else:
                tc_dict[key] = sc.UnityCurve()

        tc_master = sc.UnityCurve(lam=self.lam, lam_res=self.lam_res,
                                  min_step=self.cmds["SIM_SPEC_MIN_STEP"])
        for key in tc_keywords:
            tc_master *= tc_dict[key]

        return tc_master


    def _gen_master_psf(self, psf_type="Airy"):
        """
        Generate a Master PSF for the system. This includes the AO PSF.
        Notes: Jitter can be applied to detector array as a single PSF, and the
               ADC shift can be applied to each layer of the psf separately

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

        if self.cmds["SCOPE_PSF_FILE"] is not None and not \
                                    os.path.exists(self.cmds["SCOPE_PSF_FILE"]):
            warnings.warn("There is no PSF file under " + self.cmds["SCOPE_PSF_FILE"])


        if self.cmds["SCOPE_PSF_FILE"] is not None and \
                                os.path.exists(self.cmds["SCOPE_PSF_FILE"]):
            print("Using PSF:", self.cmds["SCOPE_PSF_FILE"])
            psf_m1 = psf.UserPSFCube(self.cmds["SCOPE_PSF_FILE"],
                                     self.lam_bin_centers)
            
            if psf_m1[0].pix_res != self.pix_res:
                psf_m1.resample(self.pix_res)

            if self.cmds["SCOPE_STREHL_RATIO"] < 1. and \
                                    "PSF_POPPY" in self.cmds["SCOPE_PSF_FILE"]:
                strehl = self.cmds["SCOPE_STREHL_RATIO"]
                gauss = psf.GaussianPSF(fwhm=0.8,
                                        size=psf_m1[0].size,
                                        pix_res=self.pix_res,
                                        undersized=True)
                psf_m1 = psf_m1 + [gauss.array * (1-strehl)]*len(psf_m1)

        else:
            print("Using PSF:", psf_type)
            ao_eff = self.cmds["SCOPE_AO_EFFECTIVENESS"]

            # Get a Diffraction limited PSF
            fwhm = (1.22*u.rad * self.lam_bin_centers * u.um / \
                                (self.cmds.diameter * u.m)).to(u.arcsec).value
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
                    psf_seeing = psf.GaussianPSF(self.lam_bin_centers,
                                                 fwhm=fwhm,
                                                 pix_res=self.pix_res)

                    psf_m1 = psf_diff.convolve(psf_seeing)
                else:
                    psf_m1 = psf_diff

        psf_master = psf_m1
        return psf_master

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
def get_filter_curve(filtername):
    """
    Return a Vis/NIR broadband filter TransmissionCurve object

    Notes
    -----
    Acceptable filter names are B V R I Y z J H K Ks
    """
    return filter_curve(filtername)


def filter_curve(filtername):
    """
    Return a Vis/NIR broadband filter TransmissionCurve object

    Notes
    -----
    Acceptable filter names are B V R I Y z J H K Ks
    """
    if filtername not in "BVRIYzJHKKs":
        raise ValueError("filter not recognised: "+filter)
    fname = os.path.join(__pkg_dir__, "data", "TC_filter_"+filtername+".dat")
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


def package_dir():
    """
    Return the folder where all the data files are stored
    """
    return __pkg_dir__


def poppy_eelt_psf_cube(lam_bin_centers, filename=None, **kwargs):
    """
    Generate a FITS file with E-ELT PSFs for a range of wavelengths by using the
    POPPY module - https://pythonhosted.org/poppy/

    Parameters:
    -----------
    lam_bin_centers: [um] the centres of each wavelength bin in micron
    filename: [None] path name to where the FITS file should be saved. If not
              given, a HDUList is returned

    Optional Parameters:
    --------------------
    pix_res: [arcsec] the angular resolution of the pixels. Default is 1 mas
    diameter_out: [m]
    diameter_in: [m]
    flatotflat: [m]
    gap: [m]
    n_spiders: [m]
    size: [int]
    oversample: [int]
    clobber: [True/False]
    """

    try:
        import poppy
    except ImportError:
        print("Please install poppy \n >>sudo pip3 install poppy")
        return

    params = {"diameter_out"  :37,
              "diameter_in"   :5.5,
              "flattoflat"    :1.45,
              "gap"           :0.004,
              "n_spiders"     :6,
              "pix_res"       :0.001,
              "size"          :255,
              "oversample"    :1,
              "clobber"       :True,
              "cpus"          :-1}
    params.update(kwargs)

    # Create the Optical Train
    rings = int(0.65 * params["diameter_out"] / params["flattoflat"])

    m1 = poppy.MultiHexagonAperture(rings=rings,
                                    flattoflat=params["flattoflat"],
                                    gap=params["gap"])
    pri = poppy.CircularAperture(radius=params["diameter_out"]/2)
    sec = poppy.SecondaryObscuration(secondary_radius=params["diameter_in"]/2,
                                     n_supports=params["n_spiders"],
                                     support_width=0.5)
    eelt = poppy.CompoundAnalyticOptic(opticslist=[m1, pri, sec], name='E-ELT')

    # Create a "detector"
    osys = poppy.OpticalSystem()
    osys.addPupil(eelt)
    osys.add_detector(pixelscale=params["pix_res"],
                      fov_arcsec=params["pix_res"] * params["size"],
                      oversample=params["oversample"])

    # list of wavelengths for which I should generate PSFs
    lam_bin_centers = np.array(lam_bin_centers)

    # Generate the PSFs in multiple threads
    try:
        import multiprocessing as mp
        if params["cpus"] < 0:
            num_cpu = max(1, mp.cpu_count() - 1)
        else:
            num_cpu = min(params["cpus"], mp.cpu_count() - 1)

        pool = mp.Pool(processes=num_cpu)
        pond = [pool.apply_async(_get_poppy_psf, (osys, lam)) \
                for lam in lam_bin_centers]
        psf_hdu = [duck.get() for duck in pond]

    except (ImportError, mp.ProcessError):  # anything else to catch? (OC)
        print("Not possible to use multiple threads")
        psf_hdu = []
        for lam in lam_bin_centers:
            psf_hdu += [_get_poppy_psf(osys, lam)]
            print(lam)

    for psfi in psf_hdu:   # was 'psf', conflict with module
        psfi.data = psfi.data.astype(np.float32)
        psfi.header["CDELT1"] = (params["pix_res"],
                                 "[arcsec] - Pixel resolution")
        psfi.header["CDELT2"] = (params["pix_res"],
                                 "[arcsec] - Pixel resolution")
        psfi.header["PSF_TYPE"] = ("POPPY", "Generated by Poppy")
        # Convert wavelength in header to [um] - originally in [m]
        psfi.header["WAVE0"] = (psfi.header["WAVE0"]*1E6,
                                "Wavelength of PSF [um]")

    if len(psf_hdu) > 1:
        psf_hdu = [psf_hdu[0]] + \
                 [fits.ImageHDU(i.data, i.header) for i in psf_hdu[1:]]

    hdulist = fits.HDUList(psf_hdu)

    if filename is None:
        return hdulist
    else:
        hdulist.writeto(filename, clobber=params["clobber"], checksum=True)


def _get_poppy_psf(osys, lam):
    """
    A self contained function for "reading out" the detector, so that we can
    use multi-threading, if available (i.e. if it isn't a windows machine)
    """
    import datetime.datetime as dt
    t = dt.now()
    print(lam, str(t.hour)+":"+str(t.minute)+":"+str(t.second))
    return osys.calcPSF(lam * 1E-6)[0]
