##############################
# A script to test progress
##############################

import sys, os

import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
from matplotlib.colors import LogNorm

from astropy.io import fits
import astropy.units as u

import PSFCube as psf
import SpectralCurve as sc
import LightObject as lo
import Detector as fpa
import utils

##############################
# Pull in User Commands
##############################

self.cmds = utils.read_config("../user_commands/master.config")
fnames = [value for value in self.cmds.values()]
for fname in fnames: 
    if os.path.exists(fname):
        self.cmds.update(utils.read_config(fname))
    else:
        warnings.warn(fname+" doesn't exist")

if user_filename is not None: 
    self.cmds.update(utils.read_config(user_filename))

    
if self.cmds["OBS_OUTPUT_DIR"] == "none":
    self.cmds["OBS_OUTPUT_DIR"] = "./"

if self.cmds["OBS_OUTPUT_NAME"] == "none":
    self.cmds["OBS_OUTPUT_NAME"] = "output.fits"

# Check if a filter curve file has been give, or a standard broadband name
if self.cmds["INST_FILTER_TC"] in ["I", "z", "Y", "J", "H", "Ks", "K"]:
    self.cmds["INST_FILTER_TC"] = "../data/TC_filter_" + \
                                    self.cmds["INST_FILTER_TC"] + ".dat"


########################################
# Set lam_bin_edges, lam_bin_centers
########################################
                                    
# if SIM_USE_FILTER_LAM is true, then use the filter curve to set the
# wavelength boundaries where the filter is < SIM_FILTER_THRESHOLD
tc_filt = sc.TransmissionCurve(self.cmds['INST_FILTER_TC'])

if self.cmds["SIM_USE_FILTER_LAM"].lower() == "yes":
    mask = np.where(tc_filt.val > self.cmds["SIM_FILTER_THRESHOLD"])[0]
    lam_min, lam_max = tc_filt.lam[mask[0]], tc_filt.lam[mask[-1]]
    self.lam_bin_edges = np.arange( lam_min, lam_max+1E-7, 
                                            self.cmds["SIM_LAM_PSF_BIN_WIDTH"])
else:
    lam_min, lam_max = self.cmds["SIM_LAM_MIN"], self.cmds["SIM_LAM_MAX"]
    self.lam_bin_edges = np.arange( tc_filt.lam[i0], tc_filt.lam[i1]+1E-7, 
                                            self.cmds["SIM_LAM_TC_BIN_WIDTH"])

self.lam_bin_centers = 0.5 * (self.lam_bin_edges[1:] + self.lam_bin_edges[:-1]) 



########################################
# Make an optical train
########################################
















                   