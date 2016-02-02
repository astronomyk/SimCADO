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
		
		self.lam_bin_centers = None
		self.lam_bin_edges   = None
		self.make_lam_bins()
		
		
		
	def update(self, filename)	
		update_config(filename, self.cmds)
		
    def __getitem__(self, kw):
        return self.cmds[kw]
	
	
	def make_lam_bins(self):
	
		if config_dict["OBS_USE_FILTER_LAM"].lower() == "yes":
			try:
				tc_filter = Throughput(self.cmds["INST_TC_FILTER"])
				tc_cutoff = self.cmds["OBS_FILTER_THRESHHOLD"]
			
				obs_lam_min = tc_filter.lam[np.where(tc_filter.val.value>tc_cutoff)[0][0]].value  - 2 * tc_filter.res.value
				obs_lam_max = tc_filter.lam[np.where(tc_filter.val.value>tc_cutoff)[0][-1]].value + 2 * tc_filter.res.value
				if self.verbose: print("Using OBS_LAM_MIN, OBS_LAM_MAX from filter curve:", obs_lam_min, obs_lam_max)
			except: 
				print "TC_FILTER: \n", tc_filter.val.value
				txt = """Filter curve contains points lower than OBS_FILTER_THRESHHOLD
				, please specify OBS_LAM_MIN and OBS_LAM_MAX in config file"""
				raise ValueError(txt)
		else:
			obs_lam_min = config_dict["OBS_LAM_MIN"]
			obs_lam_max = config_dict["OBS_LAM_MAX"]
		obs_lam_res = config_dict["OBS_LAM_BIN_WIDTH"]
		obs_pix_res = config_dict["OBS_PIXEL_SCALE"]*u.mas

		lam_bin_edges  = np.arange(obs_lam_min, obs_lam_max+obs_lam_res, obs_lam_res)*u.um
		lam_bin_centre = (lam_bin_edges[1:]+lam_bin_edges[:-1])*0.5
		if self.verbose: print "Spectral bin edges:", lam_bin_edges