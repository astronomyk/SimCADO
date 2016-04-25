"""commands.py"""
###############################################################################
# commands
#
# DESCRIPTION
#
# commands is essentially a dictionary that holds all the variables that
# the user may wish to change. It also has some set variable like 'self.pix_res'
# that can be accessed directly, instead of from the dictionary.
#
# Special parameters
# ==================
# self.lam_res
# self.lam
#
# self.lam_psf_res
# self.lam_bin_edges
# self.lam_bin_centers
# self.pix_res
# self.fpa_res
#
# self.exptime
# self.diameter
# self.area
# self.verbose
#
# TODO list
# =========
# - put in a defaults generator
# - remove the need for the user_commands directory by using internal defaults


import os, warnings
import numpy as np
try:
    import SimCADO.spectral as sc
    import SimCADO.utils as utils
except:
    import spectral as sc
    import utils as utils


class UserCommands(object):
    """
    A dictionary class that holds all the keywords from the config files plus
    a couple of most used variables, i.e. pix_res, lam_bin_edges

    Parameters
    ==========
    - user_filename: path to the user's .config file
    - master_filename: path to the master.config file
    """

    def __init__(self, user_filename=None,
                 default_filename="../user_commands/default.config"):
                 #master_filename="../user_commands/master.config"):

        #self.cmds = utils.read_config(master_filename)

        # need to generate a list, because the cmds dict will be updated
        # fnames = [value for value in self.cmds.values()]
        # for fname in fnames:
            # if os.path.exists(fname):
                # self.cmds.update(utils.read_config(fname))
            # else:
                # warnings.warn(fname+" doesn't exist.")

        self.cmds = utils.read_config(default_filename)

        if user_filename is not None:
            self.cmds.update(utils.read_config(user_filename))

        #self.cmds["CONFIG_MASTER"] = master_filename
        self.cmds["CONFIG_USER"]   = user_filename
        self.cmds["CONFIG_DEFAULT"] = default_filename

        # check the output path directory and file name
        if self.cmds["OBS_OUTPUT_DIR"] is None:
            self.cmds["OBS_OUTPUT_DIR"] = "./"

        if self.cmds["OBS_OUTPUT_NAME"] is None:
            self.cmds["OBS_OUTPUT_NAME"] = "output.fits"

        if self.cmds["SIM_PSF_OVERSAMPLE"] == "yes":
            self.cmds["PSF_MODE"] = "oversample"
        else:
            self.cmds["PSF_MODE"] = "linear_interp"

        # Check if for a filter curve file or a standard broadband name
        if self.cmds["INST_FILTER_TC"] in ["I", "z", "Y", "J", "H", "Ks", "K"]:
            self.cmds["INST_FILTER_TC"] = "../data/TC_filter_" + \
                                            self.cmds["INST_FILTER_TC"] + ".dat"

        self._update_attributes()

        if self.verbose:
            print("Read in parameters from \n"+"\n".join(fnames))


    def _update_attributes(self):
        """
        Update the UserCommand convenience attributes so that they are the same
        as the dictionary.
        """

        self.fpa_res    = self.cmds["SIM_DETECTOR_PIX_SCALE"]
        self.pix_res    = self.fpa_res / self.cmds["SIM_OVERSAMPLING"]

        # if SIM_USE_FILTER_LAM is true, then use the filter curve to set the
        # wavelength boundaries where the filter is < SIM_FILTER_THRESHOLD

        if self.cmds["SIM_USE_FILTER_LAM"].lower() == "yes":
            tc_filt = sc.TransmissionCurve(filename=self.cmds['INST_FILTER_TC'])
            mask = np.where(tc_filt.val > self.cmds["SIM_FILTER_THRESHOLD"])[0]
            i0 = np.max((mask[0] - 1, 0))
            i1 = np.min((mask[-1] + 1 , len(tc_filt.lam) - 1))
            lam_min, lam_max   = tc_filt.lam[i0], tc_filt.lam[i1]
        else:
            lam_min, lam_max = self.cmds["SIM_LAM_MIN"], self.cmds["SIM_LAM_MAX"]

        self.lam_res = self.cmds["SIM_LAM_TC_BIN_WIDTH"]
        self.lam = np.arange(lam_min, lam_max + 1E-7, self.lam_res)

        #self.lam_psf_res = self.cmds["SIM_LAM_PSF_BIN_WIDTH"]
        #self.lam_bin_edges = np.arange(lam_min,
        #                               lam_max + self.lam_psf_res + 1E-7,
        #                               self.lam_psf_res)
        # make lam_bin_edges according to how great the ADC offsets are
        self.lam_bin_edges = self._get_lam_bin_edges(lam_min, lam_max)
        self.lam_bin_centers = 0.5 * (self.lam_bin_edges[1:] + \
                                      self.lam_bin_edges[:-1])

        self.exptime    = self.cmds["OBS_EXPTIME"]
        self.diameter   = self.cmds["SCOPE_M1_DIAMETER_OUT"]
        self.area       = np.pi / 4 * (self.diameter**2 - \
                                        self.cmds["SCOPE_M1_DIAMETER_IN"]**2)

        self.cmds["SIM_N_MIRRORS"] = self.cmds["SCOPE_NUM_MIRRORS"] + \
                                     self.cmds["INST_NUM_MIRRORS"] + \
                                     self.cmds["INST_NUM_EXT_MIRRORS"]

        for key,value in zip(self.cmds.keys(), self.cmds.values()):
            if type(value) == str and value.lower() == "none":
                self.cmds[key] = None

        self.verbose = True   if self.cmds["SIM_VERBOSE"] == "yes"   else False

        self._split_categories()


    def _get_lam_bin_edges(self, lam_min, lam_max):
        """
        Generates lam_bin_edges
        """
        effectiveness = self.cmds["INST_ADC_PERFORMANCE"] / 100.
        if effectiveness == 1.:
            lam_bin_edges = np.array([lam_min, lam_max])
            return lam_bin_edges
            
        shift_threshold = self.cmds["SIM_ADC_SHIFT_THRESHOLD"]

        ## get the angle shift for each slice
        lam = np.arange(lam_min, lam_max + 1E-7, 0.001)
        angle_shift = utils.atmospheric_refraction( lam,
                                                    self.cmds["OBS_ZENITH_DIST"],
                                                    self.cmds["ATMO_TEMPERATURE"],
                                                    self.cmds["ATMO_REL_HUMIDITY"],
                                                    self.cmds["ATMO_PRESSURE"],
                                                    self.cmds["SCOPE_LATITUDE"],
                                                    self.cmds["SCOPE_ALTITUDE"])

        ## convert angle shift into number of pixels
        ## pixel shifts are defined with respect to last slice
        rel_shift = (angle_shift - angle_shift[-1]) / self.pix_res
        if np.max(np.abs(rel_shift)) > 1000:
            raise ValueError("Pixel shifts too great (>1000), check units")

        ## Rotate by the paralytic angle
        rel_shift *= (1. - effectiveness)
        int_shift = np.array(rel_shift / shift_threshold, dtype=np.int)
        idx = [np.where(int_shift == i)[0][0]
               for i in np.unique(int_shift)[::-1]]
        lam_bin_edges = np.array(lam[idx + [len(lam)-1]])

        return lam_bin_edges


    def _split_categories(self):
        self.obs    = {i:self.cmds[i] for i in self.cmds.keys() if "OBS" in i}
        self.sim    = {i:self.cmds[i] for i in self.cmds.keys() if "SIM" in i}
        self.atmo   = {i:self.cmds[i] for i in self.cmds.keys() if "ATMO" in i}
        self.scope  = {i:self.cmds[i] for i in self.cmds.keys() if "SCOPE" in i}
        self.inst   = {i:self.cmds[i] for i in self.cmds.keys() if "INST" in i}
        self.fpa    = {i:self.cmds[i] for i in self.cmds.keys() if "FPA" in i}
        self.hxrg   = {i:self.cmds[i] for i in self.cmds.keys() if "HXRG" in i}

    def update(self, x_dict):
        if isinstance(x_dict, commands):
            self.cmds.update(x_dict.cmds)
        elif isinstance(x_dict, dict):
            self.cmds.update(x_dict)
        else:
            raise ValueError("Cannot update with type: "+type(x_dict))
        self._update_attributes()

    def keys(self):
        return self.cmds.keys()

    def values(self):
        return self.cmds.values()

    def __str__(self):
        if self.cmds["CONFIG_USER"] is not None:
            return "A dictionary of commands compiled from " + \
                                                        self.cmds["CONFIG_USER"]
        else:
            return "A dictionary of default commands"

    def __iter__(self):
        return self.cmds.__iter__()

    def __getitem__(self, key):
        return self.cmds[key]

    def __setitem__(self, key, val):
        self.cmds[key] = val
        self._update_attributes()


    ### Add to the update that all the cmds.variable are updated when
    ### the dicts are updated

class bloedsinn():
    pass