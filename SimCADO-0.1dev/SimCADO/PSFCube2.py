###############################################################################
# PSFCube
#
# DESCRIPTION
#
# "If you want to bake an apple pie from scratch, 
#     first you must create the universe" - Carl Sagan
#
# === SINGLE PSFS ===
#
# We need to start by generating a single PSF in order to generate a PSFCube.
# We need to know the spatial characteristics of the PSF:
# The commonalities of all PSFs are:
#   - pix_width
#   - pix_height
#   - pix_res
#   - type
#
# The types of PSF offered: Moffat, Gaussian2D, Airy, Delta, Line, User
# For each of the PSF types we need to create a subclass of PSF. Each subclass
# takes its own list of parameters:
#   - MoffatPSF      (alpha, beta)
#   - GaussianPSF    (fwhm, eccentricity=0, angle=0)   
#   - AiryPSF        (first_zero, eccentricity=0, angle=0)
#   - DeltaPSF       (x=0, y=0)
#   - LinePSF        (x0, x1, y0, y1, angle=0)
#   - UserPSF        (filename, ext_no=0)
#
# 
# === MULTIPLE PSFS IN A CUBE ===
#
# To generate a PSF cube we need to know the spectral bins and the type of PSF. 
# The bins are defined by a central wavelength, however a cube should also
# contain the edges of each bin so that transmission and emission can be 
# re-binned properly.
#   - lam_bin_centers
#   - lam_bin_edges
#   - lam_res
#
#
# A PSF instance will have these additional arguments:
#   - array ... a 2D array to hold the PSF image
# A PSFCube instance will have these additional arguments:
#   - cube ... a (l,x,y) 3D array to hold the PSF cube
#
# As far as input goes, PSFCube should be able to accept a dictionary with the
# keywords necessary to build the cube.
#
# NOTES: 
# All wavelength values are given in [µm]
# All pixel dimensions are given in [arcsec]
# All angles are given in [deg]
#
#
# CLASSES:
#  PSF(object)
#  PSFCube(object)
#
#
# SUBCLASSES:
#  MoffatPSF(PSF)
#  GaussianPSF(PSF)
#  AiryPSF(PSF)
#  DeltaPSF(PSF)
#  LinePSF(PSF)
#  UserPSF(PSF)
#
#  PSFCubeADC(PSFCube)
#  PSFCubeFromFile(PSFCube)
#  PSFCubeAnalytic(PSFCube)
#  
#
# METHODS:
#


# There are two types of psf object here:
#     - a cube
#     - a single psf image
# The cube is essentially a list of psf images, homogenized in size
# Should we have separate classes for these?
#
# Both PSF and PSFCube can be created from a single model or through
# convolution of a list of PSF components

import numpy as np

from astropy import units as u
from astropy.convolution import (Gaussian2DKernel, AiryDisk2DKernel,
                                 Moffat2DKernel)
from astropy.convolution import convolve, convolve_fft
import utils


## These classes and functions are exported to the package
__all__ = ["PSFCube", "PSF", "MoffatPSF", "AiryPSF", "GaussianPSF", "DeltaPSF"]        


        
class PSF(object):
    """Point spread function (single layer) base class

    Needed keywords arguments:
    - size: the size length of the array in pixels
    - pix_res: the pixel scale used in the array
    """
    
    def __init__(self, size, pix_res):
    
        self.size = size
        self.pix_res = pix_res
        self.array = np.zeros((self.size, self.size))
        
        self.info = dict([])
        self.info['created'] = 'yes'
        self.info['description'] = "Point spread function (single layer)"
    
    def __repr__(self):
        return self.info['description']

    def set_array(sefl, array, threshold=1e-15):
        """
        Renormalise the array and make sure there aren't any negative values
        that will screw up the flux later on
        
        Keywords:
        - array
        - threshold: by default set to 1E-15
        """
        self.array = array
        self.array[self.array <=0] = threshold
        self.array = self.array / np.sum(self.array)
        self.size = self.array.shape[0]
        
    def convolve(self, kernel): 
        """
        Convolve the PSF with another kernel. The PSF.array keeps its shape
        - kernel is a PSF object
        - ## TODO the option to convolve with a 2D ndarray ##
        - ## TODO resize the kernel if the pixel scales are different ##
        """
        self.set_array(convolve_fft(self.array, kernel.array))

    def add():
        pass
    
    def mult():
        pass
    
    
    def resize(self, new_size):
        """Reshape the PSF

        The target shape is new_size x new_size
        """
        # make sure the new size is always an odd number
        if new_size % 2 == 0:
            new_size += 1

        new_arr = np.zeros((new_size, new_size))
        aw, bw = self.array.shape[0], new_size

        # if the new array is larger, place old array in the centre of the 
        # empty array. Otherwise, cut the centre out of the old array
        if bw > aw:
            new_arr[bw/2 - aw/2 : bw/2 + aw/2 + 1, 
                    bw/2 - aw/2 : bw/2 + aw/2 + 1] = self.array
        else:
            new_arr = self.array[aw/2 - bw/2 : aw/2 + bw/2 + 1, 
                                 aw/2 - bw/2 : aw/2 + bw/2 + 1]
        
        self.set_array(new_arr)
    
    
    
    
class DeltaPSF(PSF):
    """
    Generate a PSF with a delta function at position (x,y)
    
    Needed keywords arguments:
                
    Optional keywords
    - position: where (x,y) on the array will the delta function go,
                default is (x,y) = (0,0) and is the centre of the array
    - size: the side length of the array in pixels
    - pix_res: the pixel scale used in the array, default is 0.004 arcsec   
    """
        
    def __init__(self, **kwargs):

        if "position" not in kwargs.keys():
            self.position = kwargs["position"]  
        else: 
            self.position = (0,0)
        
        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 1
        else: 
            size = int(np.max(np.abs(position))) * 2 + 1
        
        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]  
        else: 
            pix_res = 0.004
        
        
        super(DeltaPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Delta"
        self.info['description'] = "Delta PSF, centred at (%.1f, %.1f)" \
                                    % position
        
        self.x = self.size // 2 + position[0]
        self.y = self.size // 2 + position[1]

        x2 = self.x - int(self.x)
        x1 = 1. - x2
        y2 = self.y - int(self.y)
        y1 = 1. - y2

        arr = np.zeros((self.size,self.size))
        arr[int(self.y) : int(self.y) + 2, int(self.x) : int(self.x) + 2] = \
            np.array([[x1 * y1, x2 * y1], [x1 * y2, x2 * y2]])
        self.set_array(arr)
        
 
        
class AiryPSF(PSF):
    """
    Generate a PSF for an Airy function with an equivalent FWHM
    
    Needed keywords arguments:
    - fwhm: the equivalent FWHM in [arcsec] of the PSF.
    
    Optional keywords
    - size: the side length of the array in pixels
    - pix_res: the pixel scale used in the array, default is 0.004 arcsec   
    """
    
    def __init__(self, fwhm, **kwargs):
        
        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]  
        else: 
            pix_res = 0.004

        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 1
        else: 
            size = 1
        size = int(np.max((round(8 * self.fwhm / pix_res) * 2 + 1, self.size)))

        super(AiryPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Airy"
        self.info['description'] = "Airy PSF, FWHM = %.1f arcsec" \
                                    % (self.fwhm)
                                    
        ## convert sigma (gauss) to first zero (airy)
        gauss2airy = 2.76064 
                               
        n = (self.fwhm / 2.35) / self.pix_res * gauss2airy
        self.set_array(AiryDisk2DKernel(n, x_size=self.size, y_size=self.size,
                                        mode='oversample').array)
    
class GaussianPSF(PSF):
    """
    Generate a PSF for an Gaussian function
    
    Needed keywords arguments:
    - fwhm: the FWHM in [arcsec] of the PSF.
    
    Optional keywords
    - size: the side length of the array in pixels
    - pix_res: the pixel scale used in the array, default is 0.004 arcsec   
    """
    
    def __init__(self, fwhm, **kwargs):
        
        self.fwhm = fwhm
        
        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]  
        else: 
            pix_res = 0.004

        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 1
        else: 
            size = 1
        size = int(np.max((round(8 * self.fwhm / pix_res) * 2 + 1, size)))

        super(GaussianPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Gaussian"
        self.info['description'] = "Gaussian PSF, FWHM = %.1f arcsec" \
                                    % (self.fwhm)
                                                                  
        n = (self.fwhm / 2.35) / self.pix_res
        self.set_array(Gaussian2DKernel(n, x_size=self.size, y_size=self.size,
                                        mode='oversample').array)
        
        
class MoffatPSF(PSF):
    """
    Generate a PSF for a Moffat function. Alpha is generated from the FWHM and 
    Beta = 4.765 (from Trujillo et al. 2001)
    
    Needed keywords arguments:
    - fwhm: the FWHM in [arcsec] of the PSF.
    
    Optional keywords
    - size: the side length of the array in pixels
    - pix_res: the pixel scale used in the array, default is 0.004 arcsec   
    """ 
        
     def __init__(self, fwhm, **kwargs):
        
        self.fwhm = fwhm
        
        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]  
        else: 
            pix_res = 0.004

        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 1
        else: 
            size = 1
        size = int(np.max((round(8 * self.fwhm / pix_res) * 2 + 1, size)))

        super(GaussianPSF, self).__init__(size, pix_res)

    
        beta = 4.765 ### Trujillo et al. 2001
        alpha = self.fwhm/(2 * np.sqrt(2**(1/beta) - 1))
        self.info["Type"] = "Moffat"
        self.info['description'] = "Moffat PSF, FWHM = %.1f, alpha = %.1f"\
                                       % (self.fwhm, alpha)

        if self.size > 100:
            mode = "linear_interp"
        else:
            mode = "oversample"
        self.set_array(Moffat2DKernel(alpha, beta, x_size=self.size, 
                                      y_size=self.size, mode=mode).array)
 

class CombinedPSF(PSF):


        
        
    @classmethod
    def __init__(self, psf_list):
        """Generate a master psf through convolution of a list of psfs"""
        import copy

        if not hasattr(psf_list, "__len__"):
            psf_list = [psf_list]

        self = copy.deepcopy(psf_list[0])
        self.info['description'] = "Master psf from list"
        self.info['PSF01'] = psf_list[0].info['description']
        for i, psf in enumerate(psf_list[1:]):
            self.psf = convolve_fft(self.psf, psf.psf)
            self.info['PSF%02d' % (i+2)] = psf.info['description']

        return self
        

        ## if lam does not have units, assume microns
        #if self.fwhm.value == 0:
            ## TODO: telescope parameters as variables!
            #self.fwhm = (1.22 * lam / (39.3 * u.m)).cgs * u.rad * (
                #206264806 * u.mas / u.rad)



       

                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 


class PSFCube(object):
    """Class holding wavelength dependent point spread function"""

    def __init__(self):
        self.info = dict([])
        self.info['created'] = 'yes'

    def __repr__(self):
        return self.info['description']

    @classmethod
    def gen_adc(self, lam_arr, config_dict, res=1*u.mas, nadir_angle=0):
        # TODO: Remove nadir angle (=parallactic angle?)
        # TODO: add info
        effectiveness = float(config_dict['INST_ADC_EFFICIENCY'])/100. ###
        lat = float(config_dict['SCOPE_LATITUDE'])
        alt = float(config_dict['SCOPE_ALTITUDE'])
        rel_hum = float(config_dict['ATMO_REL_HUMIDITY'])

        ## CHECK: better called OBS_ZENITH_DIST
        z0 = float(config_dict['OBS_ZENITH_DIST'])
        temp = float(config_dict['ATMO_TEMPERATURE'])
        pres = float(config_dict['ATMO_PRESSURE'])

        lam_arr = unify(lam_arr, u.um)
        nadir_angle = unify(nadir_angle, u.deg)
        res = unify(res, u.mas)

        ## get the angle shift for each slice
        angle_shift = [atmospheric_refraction(lam, z0, temp, rel_hum, pres,
                                              lat, alt)
                       for lam in lam_arr.value]
        angle_shift = angle_shift * u.arcsec

        ## convert angle shift into number of pixels
        ## pixel shifts are defined with respect to last slice
        pixel_shift = (angle_shift - angle_shift[-1]).to(u.mas) / res

        ## Rotate by the parallactic angle
        x = -pixel_shift * np.sin(nadir_angle.to(u.rad)) * (1. - effectiveness)
        y = -pixel_shift * np.cos(nadir_angle.to(u.rad)) * (1. - effectiveness)

        adc_cube = PSFCube.gen_cube(lam_arr, kernel='adc', position=(x, y),
                                    res=res)
        return adc_cube

    @classmethod
    def gen_cube(self, lam_arr, fwhm=0, res=1.*u.mas, kernel='gauss',
                 position=(0,0), padding=5, min_size=71, size=None):
        """Generate a cube full of psf_gen objects

        An astropy unit array is required for the central wavelength
        of each slice.
        The FWHM can be specified as a blanket value, or an array of
        FWHMs for each slice
        """

        ## Get everything in mas for the spatial values and um for the
        ## spectral values
        self = PSFCube()

        if 'airy' in kernel:
            self.info['description'] = "PSF cube, Airy"
        elif 'gauss' in kernel:
            self.info['description'] = "PSF cube, Gauss"
        else:
            self.info['description'] = "PSF cube, delta peak"
            
        self.lam_arr = utils.unify(lam_arr, u.um)

        # if fwhm is scalar, generate a constant array
        self.fwhm = utils.unify(fwhm, u.mas, len(lam_arr))
        self.res = utils.unify(res, u.mas)
        self.kernel = kernel
        self.padding = padding

        ## Create a cube along the spectral dimension. Call psf_gen for each
        ## lambda value and corresponding fwhm
        self.cube = []

        ## if only one position, generate constant arrays for x and y
        x, y = position[0], position[1]
        if not hasattr(x, "__len__"):
            x = x * np.ones(len(lam_arr))
        if not hasattr(y, "__len__"):
            y = y * np.ones(len(lam_arr))
            
        for i in range(len(self.lam_arr)):
            self.cube += [PSF.gen_analytic(self.lam_arr[i], fwhm=self.fwhm[i],
                                           res=self.res, kernel=self.kernel,
                                           position=(x[i], y[i]),
                                           padding=padding,
                                           min_size=min_size)]

        self.fwhm = [i.fwhm for i in self.cube]
        self.size_orig = [i.psf.shape[0] for i in self.cube]

        ## CHECK: Do we have to resize the slices to a common shape?
        ## Comment from Kieran's original code
        # resample the slices so that they are all the same size
        # self.resample(np.max(self.size_orig))
        
    # !!!!!!!!!!! It seems that we don't need the resampling,
        # so long as all arrays are odd numbers in length !!!!!!!!!!!!!!

        ## TODO: Add some more cube info

        return self
        

    @classmethod
    def gen_from_list(self, psf_list):
        """Generate a master psf cube through convolution of a list of psfs"""
        import copy

        if not hasattr(psf_list, "__len__"):
            self.psf_list = [psf_list]

        self = PSFCube()
        
        self.info['description'] = "Master psf cube from list"
        for i in range(len(psf_list)):
            self.info['PSF%02d' % (i+1)] = psf_list[i].info['description']
        
        ## Check that the wavelengths are equal
        lam_list = [psf.lam_arr for psf in psf_list]
        if not all([all(lam == lam_list[0]) for lam in lam_list[1:]]):
            raise ValueError("Wavelength arrays of psf cubes are not equal")
        self.lam = lam_list[0]

        ## Check that the resolutions are equal
        res_list = [psf.res for psf in psf_list]
        if not all([res == res_list[0] for res in res_list[1:]]):
            raise ValueError("Resolutions of psf cubes are not equal")
        self.res = res_list[0]

        ## Combine each layer
        # CHECK: What's this padding for? size is not defined
        ##padding = np.sum(np.max([psf.size for psf in psf_list], axis=1))

        self.psf = PSFCube.gen_cube(self.lam, res=self.res, kernel='dot')
        #, min_size=min_size, padding=padding)

        self.cube = range(len(self.lam))
        for i in range(len(self.lam)):
            self.cube[i] = PSF.gen_from_list([psf.cube[i] for psf in psf_list])
        
        return self