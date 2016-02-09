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

from copy import deepcopy
from astropy import units as u
from astropy.io import fits, ascii
import numpy as np
import warnings

__all__ = ["TransmissionCurve", "EmissionCurve"] 


class TransmissionCurve(object):
    def __init__(self, **kwargs):
        """
        Very basic class to either read in a text file for a transmission curve
        or take two vectors to make a transmission curve
        
        List of kwargs:
        lam: [µm] 1D numpy array of length n
        val: 1D numpy array of length n
        res: [µm] float with the desired spectral resolution
        filename: string with the path to the transmission curve file where
                  the first column is wavelength in [µm] and the second is the
                  transmission coefficient between [0,1]
        """
        self.params = { "lam_unit"  :u.um,
                        "val_unit"  :None,
                        "filename"  :None,
                        "lam_res"   :0.001,
                        "Type"      :"Transmission",
                        "min_bin_width" :1E-5
                       }
                       
        self.params.update(kwargs)
        
        self.info = dict([])
        self.info["Type"] = self.params["Type"]
        
        self.lam_orig, self.val_orig = self.get_lam_val()
        self.lam_orig *= (1*self.params["lam_unit"]).to(u.um)
        
        self.resample(self.params["lam_res"])

        
    def __repr__(self):
        return "Ich bin eine SpectralCurve:\n"+str(self.info)

        
    def get_lam_val(self):
        """
        Get the wavelength and value vectors from the input parameters
        """
    
        if "lam" in self.params.keys() and "val" in self.params.keys():
            lam = self.params["lam"]
            val = self.params["val"]
            
        # test if it is a skycalc file
        elif "filename" in self.params.keys():
            filename = self.params["filename"]
            if ".fits" in filename:
                hdr = fits.getheader(filename)
                if any(["SKYCALC" in hdr[i] for i in range(len(hdr)) \
                                                    if type(hdr[i]) == str]):
                    if self.params["Type"] == "Emission":
                        lam = fits.getdata(filename)["lam"]        
                        val = fits.getdata(filename)["flux"]
                    else:
                        lam = fits.getdata(filename)["lam"]
                        val = fits.getdata(filename)["trans"]
                else:
                    data = fits.getdata("../data/skytable.fits")
                    lam = data[data.columns[0].name]
                    val = data[data.columns[1].name]
            else:
                data = ascii.read(self.params["filename"])
                lam = data[data.colnames[0]]
                val = data[data.colnames[1]]
        else: 
            raise ValueError("Please pass either filename or lam/val keywords")
        
        return lam, val
        
    def resample(self, bins, action="average", use_edges=False):
        """
        Resamples both the wavelength and value vectors to an even grid. 
        In order to avoid losing spectral information, the TransmissionCurve
        resamples down to a resolution of 'min_bin_width' (default: 0.01nm) 
        before resampling again up to the given sampling vector defined by
        'bins'.
        
        Keywords:
        - bins: [µm]: float - taken to mean the width of bins on an even grid
                      array - the centres of the spectral bins
        
        Optional keywords:
        - action: ['average','sum'] How to rebin the spectral curve. If 'sum', 
                  then the curve is normalised against the integrated value of  
                  the original curve. If 'average', the average value per bin
                  becomes the value for each bin.
        - use_edges: [False, True] True if the array passed in 'bins' describes 
                     the edges of the wavelength bins. 
        
        """
        #####################################################
        # Work out the irregular grid problem while summing #
        #####################################################
        
        min_step = self.params["min_bin_width"]
        
        tmp_x = np.arange(self.lam_orig[0], self.lam_orig[-1], min_step)
        tmp_y = np.interp(tmp_x, self.lam_orig, self.val_orig)
        
        # The summing issue - assuming we want to integrate along the curve,
        # i.e. count all the photons in a new set of bins, we need to integrate
        # along the well-sampled (1E-5µm) curve. However the above line of code
        # using np.interp doesn't take into account the new bin width when 
        # resampling down to 1E-5µm. I account for this summing up all the
        # photons in the original data set and normalising the new 1E-5 bin 
        # data set to have the same amount.
        if action == "sum": tmp_y *= (np.sum(self.val_orig) / np.sum(tmp_y))
        
        # if bins is a single number, use it as the bin width
        # else as the bin centres
        if not hasattr(bins, "__len__"): 
            lam_tmp = np.arange(self.lam_orig[0], self.lam_orig[-1], bins)
        else: 
            lam_tmp = bins
        lam_res = bins[1] - bins[0]
        
        # define the edges and centres of each wavelength bin
        if use_edges:
            lam_bin_edges = bins
            lam_bin_centres = 0.5 * (bins[1:] + bins[:-1])
        else:
            lam_bin_edges = np.append(bins - 0.5*lam_res, bins[-1] + 0.5*lam_res)
            lam_bin_centers = bins
       
        # here is the assumption of a regular grid - see res_tmp
        val_tmp = np.zeros((len(lam_bin_centers)))
        
        for i in range(len(lam_bin_centers)):
            
            mask_i = np.where((tmp_x > lam_bin_edges[i]) * 
                              (tmp_x < lam_bin_edges[i+1]))[0]     
            
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
    
    
    def __getitem__(self, i):
        return self.val[i], self.lam[i]
        
    def __array__(self):
        return self.val
    
    def __mul__(self, tc):
        """ 
        Product of a TransmissionCurve with a scalar or another Curve
        If tc is a TransmissionCurve and does not have the same lam, it is 
        resampled first.
        """
        tcnew = deepcopy(self)
        
        if not hasattr(tc, "val"):
            tcnew.val *= tc
        else:
            ### TODO: This comparison needs work
            if not np.all(self.lam == tc.lam):
                tc.resample(self.lam)
            tcnew.val *= tc.val

        return tcnew

    def __add__(self, tc):
        """ 
        Addition of a TransmissionCurve with a scalar or another Curve.
        If tc is a TransmissionCurve and does not have the same lam, it is 
        resampled first.
        """
        tcnew = deepcopy(self)

        if not hasattr(tc, "val"):
            tcnew.val += tc
        else:
            ### TODO: This comparison needs work
            if not np.all(self.lam == tc.lam):
                tc.resample(self.lam)
            tcnew.val += tc.val

        return tcnew 
        
    def __sub__(self, tc):
        return  self.__add__(-1 * tc)
        
    def __rmul__(self, x):
        return self.__mul__(x)
                
    def __radd__(self, x):
        return self.__add__(x)
    
    def __rsub__(self, x):
        self.__mul__(-1)
        return self.__add__(x)         

    def __imul__(self, x):
        return self.__mul__(x)
                
    def __iadd__(self, x):
        return self.__add__(x)
    
    def __isub__(self, x):
        return self.__sub__(x)         
        
        
        
class EmissionCurve(TransmissionCurve):
    def __init__(self, **kwargs):
        """
        List of kwargs:
        - lam: 1D numpy array of length n in [µm]
        - val: 1D numpy array of length n in []
        - res: float with the desired spectral resolution in [µm]
        - filename: string with the path to the transmission curve file where
                  the first column is wavelength in [µm] and the second is the
                  transmission coefficient between [0,1]
        
        - pix_res: [arcsec] float of int for the field of view for each pixel
        - area: [m2] float or int for the collecting area of M1
        - exptime: [s] float or int for the integration time for an exposure
        - units: string or astropy.units for calculating the number of photons 
               per voxel
        """
        default_params = {  "pix_res" :0.004,
                            "area"    :978,
                            "exptime" :1,
                            "val_unit":"ph/(s m2 micron arcsec2)"
                          }

        if "val_unit" not in kwargs.keys():
            warnings.warn("No val_unit specified in EmissionCurve. Assuming ph/(s m2 micron arcsec2)")
                          
        super(EmissionCurve, self).__init__(Type = "Emission", **kwargs)
        self.params.update(default_params)
        self.convert_to_photons()
        
    def resample(self, bins, action="sum"):
        super(EmissionCurve, self).resample(bins, action)

    def convert_to_photons(self):
        """Do the conversion to photons/voxel by using the val_unit, lam, area
        and exptime keywords. If not given, make some assumptions.
        """
        self.params["val_unit"] = u.Unit(self.params["val_unit"])
        bases  = self.params["val_unit"].bases
        
        factor = 1. * self.params["val_unit"]

        if u.s      in bases: factor *= self.params["exptime"] 
        if u.m      in bases: factor *= self.params["area"]
        if u.arcsec in bases: factor *= self.params["pix_res"]**2
        if u.micron in bases: factor *= self.params["lam_res"]

        self.val *= factor
