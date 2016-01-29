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
# TransmissionCurve rebin would involve averaging the "val" values, while a
# SpectralCurve rebin would involve summing up the "val" values.
# 
# SpectralCurve also needs the file path of the data:
#  - filename
#
# A TransmissionCurve doesn't need anything else on top of the SpectralCurve,
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
from astropy.io import fits
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
##    @classmethod
##    def from_ascii(self, fname, lamcol=0, valcol=1):
##        '''Read a spectral from an ascii file.'''
##
##        d = np.loadtxt(fname, usecols=(lamcol, valcol))
##        speccurve = SpectralCurve(d[:, 0], d[:, 1])
##        return speccurve
        

class Emission(SpectralCurve):
    """Class holding spectral emission curve"""
    def __init__(self, lam, val, Type=None):
        self = SpectralCurve(lam, val, Type=Type)

    def __repr__(self):
        return "Ich bin eine Emission"

        

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
