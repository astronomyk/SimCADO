###############################################################################
# UserCommands
#
# DESCRIPTION
#
# UserCommands is essentially a dictionary that holds all the variables that
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
# self.fpa_pix_res
#
# self.exptime
# self.diameter
# self.area 
# self.verbose


import os, warnings
import numpy as np
import SpectralCurve as sc
import utils

class UserCommands(object):
    """
    A dictionary class that holds all the keywords from the config files plus
    a couple of most used variables, i.e. pix_res, lam_bin_edges
    
    Parameters
    ==========
    - user_filename: path to the user's .config file 
    - master_filename: path to the master.config file
    """
    
    def __init__(self, user_filename=None, master_filename="../user_commands/master.config"):

        self.cmds = utils.read_config(master_filename)
        
        # need to generate a list, because the cmds dict will be updated
        fnames = [value for value in self.cmds.values()]
        for fname in fnames: 
            if os.path.exists(fname):
                self.cmds.update(utils.read_config(fname))
            else:
                warnings.warn(fname+" doesn't exist.")

        if user_filename is not None: 
            self.cmds.update(utils.read_config(user_filename))   
        
        self.cmds["CONFIG_MASTER"] = master_filename
        self.cmds["CONFIG_USER"]   = user_filename
        
        # check the output path directory and file name
        if self.cmds["OBS_OUTPUT_DIR"] == "none":
            self.cmds["OBS_OUTPUT_DIR"] = "./"
   
        if self.cmds["OBS_OUTPUT_NAME"] == "none":
            self.cmds["OBS_OUTPUT_NAME"] = "output.fits"
   
        if self.cmds["SIM_PSF_OVERSAMPLE"] == "yes":
            self.cmds["PSF_MODE"] = "oversample"
        else:
            self.cmds["PSF_MODE"] = "linear_interp"
   
        # Check if for a filter curve file or a standard broadband name
        if self.cmds["INST_FILTER_TC"] in ["I", "z", "Y", "J", "H", "Ks", "K"]:
            self.cmds["INST_FILTER_TC"] = "../data/TC_filter_" + \
                                            self.cmds["INST_FILTER_TC"] + ".dat"

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
        
        self.lam_psf_res = self.cmds["SIM_LAM_PSF_BIN_WIDTH"]
        self.lam_bin_edges = np.arange(lam_min, 
                                       lam_max + self.lam_psf_res, 
                                       self.lam_psf_res)
        self.lam_bin_centers = 0.5 * (self.lam_bin_edges[1:] + \
                                      self.lam_bin_edges[:-1])
        
        self.pix_res     = self.cmds["SIM_INTERNAL_PIX_SCALE"]
        self.fpa_pix_res = self.cmds["SIM_DETECTOR_PIX_SCALE"]
        
        self.exptime     = self.cmds["OBS_EXPTIME"]
        self.diameter    = self.cmds["SCOPE_M1_DIAMETER_OUT"]
        self.area        = np.pi / 4 * (self.diameter**2 - \
                                        self.cmds["SCOPE_M1_DIAMETER_IN"]**2)
        
        self.verbose = True     if self.cmds["VERBOSE"] == "yes"    else False
        
        if self.verbose:
            print("Read in parameters from ")
    
    def __repr__(self):
        return "A dictionary of commands compiled from "+self.cmds["CONFIG_MASTER"]
    
    def __getitem__(self, kw):
        return self.cmds[kw]
        
    def __setitem__(self, kw, val):
        self.cmds[kw] = val
        
    def keys(self):
        return self.cmds.keys()
        
    def values(self):
        return self.cmds.values()