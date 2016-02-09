###############################################################################
# SpectralCurve
#
# DESCRIPTION
#
# "If you pay peanuts, you get monkeys"
#
# A SpectralCurve is the base class for either a transmission curve or an 
# emission curve. The main attributes are 2 equal length arrays holding the 
# centres of each wavelength bin and the corresponding value - an energy or a 
# transmission factor [0-1]
#  - lam
#  - val 
#
# SpectralCurve should be overloaded on the + and * operators. Although a rebin
# method would be good, this is unique to the subclasses. E.g. a 
# Throughput curve rebin would involve averaging the "val" values, while a
# SpectralCurve rebin would involve summing up the "val" values.
# 
# SpectralCurve also needs the file path of the data:
#  - filename
#
# A Throughput curve doesn't need anything else on top of the SpectralCurve,
# however the EmissionCurve must know which units are being used so that it
# can immediately convert the energy into photons. In order to do this the 
# EmissionCurve needs the following extra info from the UserCommands dictionary
#  - spatial area [m2]
#  - angular area [arcsec2]
#  - integration time [s]
#
# As stated above each subclass should have its own rebin(lam_res) method 
#
# Notes:
# All wavelength values are in [µm]
# All other values are either transmission [0-1] or number of photons [>=0]
#
# Classes:
#  SpectralCurve(object) 
#  - from_file(filename)
#  - from_list([ThroughputCurve])
#  - from_skycalc(filename)
#
# Subclasses:
# Emission(SpectralCurve)
# - rebin(lam)
# Throughput(SpectralCurve)
# - rebin(lam)
#
# Methods:
#	
#
#


from astropy import units as u
from astropy.io import fits, ascii
import numpy as np

__all__ = ["SpectralCurve", "Emission", "Throughput"] 

class SpectralCurve(object):
	def __init__(self, lam, val, lam_unit=u.um, val_unit=1., Type=None):
		# Can be created from a single data file or from a list of Throughput
		# instances
		self.info = dict([])
		self.info["Type"] = Type
		self.lam = lam
		self.lam_unit = lam_unit
		self.val = val
		self.val_unit=val_unit

	def __repr__(self):
		return "Ich bin eine SpectralCurve:\n"+str(self.info)

	def __mul__(self, tc):
		""" Product of a spectral curve with a scalar or another throughput

If tc is a throughput and does not have the same lam, it is resampled first."""
		tcnew = self
		
		if not hasattr(tc, "val"):
			tcnew.val *= tc
		else:
			### TODO: This comparison needs work
			if not np.all(self.lam == tc.lam):
				tc.resample(self.lam)
			tcnew.val *= tc.val

		return tcnew
			

	def __imul__(self, tc):
		"""x.__imul__(y) <==> x *= y"""
		self = self * tc
		return self

## CHECK: Needed?
##	  @classmethod
##	  def from_ascii(self, fname, lamcol=0, valcol=1):
##		  '''Read a spectral from an ascii file.'''
##
##		  d = np.loadtxt(fname, usecols=(lamcol, valcol))
##		  speccurve = SpectralCurve(d[:, 0], d[:, 1])
##		  return speccurve
		

class Emission(SpectralCurve):
	"""Class holding spectral emission curve"""
	def __init__(self, lam, val, Type=None):
		self = SpectralCurve(lam, val, Type=Type)

	def __repr__(self):
		return "Ich bin eine EmissionCurve"

			   

class Throughput(SpectralCurve):
	"""Class holding information about telescope and instrument throughput

	The class can be instantiated in various ways:
	- tc = Throughput(lam, trans)
	- tc = Throughput.from_skycalc(filename, res)
	"""

	def __init__(self, lam, val, Type=None):
		# Can be created from a single data file or from a list of Throughput
		# instances
		super(Throughput, self).__init__(lam, val, Type)

	def __repr__(self):
		result = "Throughput\nType: " + self.info['Type']
		return result
		
	@classmethod
	def from_skycalc(self, fname, res=0.001*u.um):

		lam = fits.getdata(fname)["lam"]
		val = fits.getdata(fname)["trans"]
		self = Throughput(lam, val)
		self.info = dict([])
		self.info["Type"] = "Atmospheric transmission"
		self.info["Filename"] = fname
		self.info["Resolution"] = res

		return self

	@classmethod
	def from_ascii(self, fname, lam_col=0, val_col=1, Type=None):
		"""Read a transmission curve from an ascii file

		The file is assumed to have no header.
		
		Parameters:
		==========
		fname - file name
		lam_col - column with lambda values [default 0]
		val_col - column with transmission values [default 1]
		Type - string describing the file, (available through 
			   the info dictionary)
		"""
		## TODO: Generalize. Can this be a method of the superclass? 
		data = np.loadtxt(fname)
		tc = Throughput(data[:,lam_col], data[:,val_col])
		tc.info['Type'] = Type
		return tc
		
	@classmethod
	def from_throughput(self, tclist):
		"""Create a throughput by multiplying a list of throughputs"""

		import copy
		
		if not hasattr(tclist, "__len__"):
			tclist = [tclist]

		self = copy.deepcopy(tclist[0])
		self.__class__ = tclist[0].__class__
		
		## TODO: this works best if the reference tc has the highest
		## resolution. Can we identify this automatically?
		for tc in tclist[1:]:
			self *= tc

		self.info['Type'] = "Master throughput"
		for i in range(len(tclist)):
			self.info['TC%02d' % (i+1)] = tclist[i].info['Type']
		return self

	def resample(self, lam):
		"""Resample to a new lambda vector by interpolating"""
		
		# Save the original attributes
		self.lam_orig = self.lam
		self.val_orig = self.val

		# If lam is a scalar, treat it as a resolution, and build
		# lam vector to cover original
		if not hasattr(lam, "__len__"):
			if type(lam) != u.quantity.Quantity:
				lam *= u.um
			lam = lam.to(self.lam_orig.unit)
			lam = np.arange(self.lam_orig.value[0],
							self.lam_orig.value[-1],
							self.lam.value) * lam.unit

		self.lam = lam
		if hasattr(lam, "value"):
			self.val = np.interp(lam.value, self.lam_orig.value, self.val_orig)
		else:
			self.val = np.interp(lam, self.lam_orig, self.val_orig)

			
			
			
			
			
# Until I find out what's going on here, I'll use a simple TransmissionCurve

class TransmissionCurve(object):
	def __init__(self, **kwargs):
		"""
		Very basic class to either read in a text file for a transmission curve
		or take two vectors to make a transmission curve
		
		List of kwargs:
		lam: 1D numpy array of length n in [µm]
		val: 1D numpy array of length n in []
		res: float with the desired spectral resolution in [µm]
		filename: string with the path to the transmission curve file where
		          the first column is wavelength in [µm] and the second is the
				  transmission coefficient between [0,1]
		"""
		
		if "lam" in kwargs.keys() and "val" in kwargs.keys():
			self.lam_orig = kwargs["lam"]
			self.lam_orig = kwargs["val"]
		if "filename" in kwargs.keys():
			data = ascii.read(kwargs["filename"])
			self.lam_orig = data[data.colnames[0]]
			self.val_orig = data[data.colnames[1]]
		if "res" in kwargs.keys():
			self.res = kwargs["res"]
		else:
			self.res = 0.001

		self.resample(self.res)
		
		
	def resample(self, bins, action="average"):
		"""
		Resamples both the wavelength and value vectors to an even grid. 
		In order to avoid losing spectral information, the TransmissionCurve
		
		"""

		#####################################################
		# Work out the irregular grid problem while summing #
		#####################################################
		
		min_step = 1E-5
		
		
		tmp_x = np.arange(self.lam_orig[0], self.lam_orig[-1], min_step)
		tmp_y = np.interp(tmp_x, self.lam_orig, self.val_orig)
		
		# The summing issue - assuming we want to integrate along the curve,
		# i.e. count all the photons in a new set of bins, we need to integrate
		# along the well-sampled (1E-5µm) curve. However the above line of code
		# using np.interp doesn't tak into account the new bin width when 
		# resampling down to 1E-5µm. I account for this summing up all the
		# photons in the original data set and normalising the new 1E-5 bin 
		# data set to have the same amount.
		if action == "sum":	tmp_y *= (np.sum(self.val_orig) / np.sum(tmp_y))
		
		
		# if bins is a single number, use it as the bin width
		# else as the bin centres

		if not hasattr(bins, "__len__"): 
			lam_tmp = np.arange(self.lam_orig[0], self.lam_orig[-1], bins)
		else: 
			lam_tmp = bins
		
		# here is the assumption of a regular grid - see res_tmp
		res_tmp = lam_tmp[1] - lam_tmp[0]
		val_tmp = np.zeros((len(lam_tmp)))
				
		for i in range(len(lam_tmp)):
			mask_i = np.where((tmp_x > lam_tmp[i] - res_tmp/2.) * 
					          (tmp_x < lam_tmp[i] + res_tmp/2.))[0]		
			
			if np.sum(mask_i) > 0 and action == "average":	
				val_tmp[i] = np.average(tmp_y[mask_i[0]:mask_i[-1]])
			
			elif np.sum(mask_i) > 0 and action == "sum":	
				# FIXED. THE SUMMING ISSUE. TEST IT         #
				# Tested - the errors are on the 0.1% level #
				val_tmp[i] = np.trapz(tmp_y[mask_i[0]:mask_i[-1]])
			
			else: 
				val_tmp[i] = 0
			
		self.lam = lam_tmp
		self.val = val_tmp
		
		
class EmissionCurve(TransmissionCurve):
	def __init__(self, **kwargs):
		"""
		List of kwargs:
		lam: 1D numpy array of length n in [µm]
		val: 1D numpy array of length n in []
		res: float with the desired spectral resolution in [µm]
		filename: string with the path to the transmission curve file where
		          the first column is wavelength in [µm] and the second is the
				  transmission coefficient between [0,1]
		
		pix_res: float of int in [arcsec] for the field of view for each pixel
		area: float or int in [m2] for the collecting area of M1
		exptime: float or int in [s] for the integration time for an exposure
		units: string or astropy.units for calculating the number of photons 
		       per voxel
		"""
		
		if "pix_res" in kwargs.keys(): self.pix_res = kwargs["pix_res"]
		#if "lam_res" in kwargs.keys(): self.lam_res = kwargs["lam_res"]
		if "area" in kwargs.keys():    self.area    = kwargs["area"]
		if "exptime" in kwargs.keys(): self.exptime = kwargs["exptime"]
		if "units" in kwargs.keys():   self.units   = kwargs["units"]
		###### Add support for both astropy units and a string with units #####
		
		super(EmissionCurve, self).__init__(**kwargs)
		
	def resample(self, bins, action="sum"):
		super(EmissionCurve, self).resample(bins, action)

	def convert_to_photons(self):
		"""Do the conversion to photons/voxel by using the units, lam, area
		and exptime keywords. If not given, make some assumptions.
		"""
		pass













		