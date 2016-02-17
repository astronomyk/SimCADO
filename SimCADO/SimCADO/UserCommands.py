###############################################################################
# UserCommands
#
# DESCRIPTION
#
# UserCommands is essentially a dictionary that holds all the variables that
# the user may wish to change, plus some methods to get certain subclasses
# out, i.e. anything that starts with "ATMO"
#
#
# Classes:
#   
#
# Methods:
#   
#
#

from utils import *
from SpectralCurve import Throughput

class UserCommands(object):

	def __init__(self, filename)

        self.filename = filename
		self.cmds = read_config(filename)
		self.verbose = True if self.cmds["VERBOSE"].lower() == "yes" else False
		
		self.cmds["LAM_BIN_CENTERS"] = None
		self.cmds["LAM_BIN_EDGES"]   = None
		self.cmds["LAM_BIN_MIN"]     = None
		self.cmds["LAM_BIN_MAX"]     = None
		self.cmds["LAM_BIN_RES"]     = None
				
		self.cmds["PIX_RES"]		 = None
		self.cmds["PIX_SHAPE"]       = None
		
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
                                              
    
    
    
    
    def update(self, filename)	
		update_config(filename, self.cmds)
	
    def __getitem__(self, kw):
        return self.cmds[kw]
                                              