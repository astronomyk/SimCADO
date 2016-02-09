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
		
		
	def update(self, filename)	
		update_config(filename, self.cmds)
		
    def __getitem__(self, kw):
        return self.cmds[kw]