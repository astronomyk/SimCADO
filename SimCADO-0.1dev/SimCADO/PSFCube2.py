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
import scipy.ndimage.interpolation as spi

from astropy.io import fits
from astropy import units as u
from astropy.convolution import (Gaussian2DKernel, AiryDisk2DKernel,
                                 Moffat2DKernel)
from astropy.convolution import convolve, convolve_fft
import utils

import warnings

## These classes and functions are exported to the package
__all__ = ["PSFCube", "PSF", "MoffatPSF", "AiryPSF", "GaussianPSF", "DeltaPSF"]        


###############################################################################
#                            PSF and PSF subclasses                           #
###############################################################################

# (Sub)Class    Needed              Optional
# PSF           size, pix_res
# DeltaPSF                          pix_res=0.004, size=f(position), position=(0,0)
# AiryPSF       fwhm                pix_res=0.004, size=f(fwhm)
# GaussianPSF   fwhm                pix_res=0.004, size=f(fwhm)
# MoffatPSF     fwhm                pix_res=0.004, size=f(fwhm)
# CombinedPSF   psf_list            size=f(psf_list)
# UserPSF       filename,           pix_res=0.004, size=f(filename), fits_ext=0


        
class PSF(object):
    """Point spread function (single layer) base class

    Needed keywords arguments:
    - size: [int] the size length of the array in pixels
    - pix_res: [arcsec] the pixel scale used in the array
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

    def set_array(self, array, threshold=1e-15):
        """
        Renormalise the array and make sure there aren't any negative values
        that will screw up the flux later on
        
        Keywords:
        - array: [ndarray] the array representing the PSF
        - threshold: by default set to 1E-15
        """
        self.array = array
        self.array[self.array <=0] = threshold
        self.array = self.array / np.sum(self.array)
        self.size = self.array.shape[0]
        
    def resize(self, new_size):
        """Resize the PSF. The target shape is (new_size, new_size).

        Keywords:
        - new_size: [int] the new size of the PSF array in pixels
        """
        
         # make sure the new size is always an odd number
        if new_size % 2 == 0:
            new_size += 1

        arr_tmp = np.zeros((new_size, new_size)) 
        arr_tmp[new_size//2, new_size//2] = 1
        self.set_array(convolve_fft(arr_tmp, self.array))

    def resample(self, new_pix_res):
        """
        Resample the PSF array onto a new grid - Not perfect, but conserves flux
        Example: new_PSF = old_PSF.resample(new_pix_res)
        
        Keywords:
        - new_pix_res: [arcsec] the pixel resolution of the returned array
        """
        scale_factor = self.pix_res / np.float(new_pix_res) 
        new_arr = spi.zoom(self.array, scale_factor, order=1)
        new_arr *= np.sum(self.array)/np.sum(new_arr) 
        
        ############################################################
        # Not happy with the way the returned type is not the same #
        # as the original type. The new object is a plain PSF      #
        ############################################################
        new_psf = PSF(size = new_arr.shape[0], pix_res = new_pix_res)
        new_psf.set_array(new_arr)
        return new_psf
        
    def convolve(self, kernel):
        """
        Convolve the PSF with another kernel. The PSF keeps its shape
        
        Keywords:
        - kernel: [PSF/ndarray] either a numpy.ndarray or a PSF (sub)class
        """       
        if issubclass(type(kernel), PSF):
            self.set_array(convolve_fft(self.array, kernel.array))
        else:
            self.set_array(convolve_fft(self.array, kernel))
            
    def __array__(self):
        return self.array
        
    def __mul__(self, x):
        self.array *= x
        return self
        
    def __add__(self, x):
        self.array += x
        return self
        
    def __sub__(self, x):
        self.array -= x
        return self
    
    def __rmul__(self, x):
        self.__mul__(x)
                
    def __radd__(self, x):
        self.__add__(x)
    
    def __rsub__(self, x):
        self.__sub__(x)        
    
    
class DeltaPSF(PSF):
    """
    Generate a PSF with a delta function at position (x,y)
    
    Needed keywords arguments:
                
    Optional keywords
    - position: [(x,y,)] where (x,y) on the array will the delta function go,
                default is (x,y) = (0,0) and is the centre of the array
    - size: [int] the side length of the array in pixels
    - pix_res: [arcsec] the pixel scale used in the array, default is 0.004 
    """
        
    def __init__(self, **kwargs):

        if "position" in kwargs.keys():
            self.position = kwargs["position"]  
        else: 
            self.position = (0,0)
        
        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 3
        else: 
            size = int(np.max(np.abs(self.position))) * 2 + 3

        if not np.max(self.position) < size:
            raise ValueError("positions are outside array borders:")
        
        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]  
        else: 
            pix_res = 0.004
        
        super(DeltaPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Delta"
        self.info['description'] = "Delta PSF, centred at (%.1f, %.1f)" \
                                    % self.position
        
        self.x = self.size // 2 + self.position[0]
        self.y = self.size // 2 + self.position[1]

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
    - fwhm: [arcsec] the equivalent FWHM of the Airy disk core.
    
    Optional keywords
    - size: [int] the side length of the array in pixels
    - pix_res: [arcsec] the pixel scale used in the array, default is 0.004 
    """
    
    def __init__(self, fwhm, **kwargs):
               
        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]  
        else: 
            pix_res = 0.004

        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 1
        else: 
            size = 255           # min_size
        
        self.fwhm = fwhm
        size = int(np.max((round(8 * self.fwhm / pix_res) * 2 + 1, size)))
        
        if size > 511: 
            size = 511
            print("FWHM [arcsec]:", fwhm, "- pixel res [arcsec]:", pix_res)
            print("Array size:", size,"x",size, "- PSF FoV:", size * pix_res)
            warnings.warn("PSF dimensions too large")
        
        super(AiryPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Airy"
        self.info['description'] = "Airy PSF, FWHM = %.1f mas" \
                                    % (self.fwhm * 1E3)
                                    
        ## convert sigma (gauss) to first zero (airy)
        gauss2airy = 2.76064 
        
        n = (self.fwhm / 2.35) / self.pix_res * gauss2airy
        self.set_array(AiryDisk2DKernel(n, x_size=self.size, y_size=self.size, \
                                        mode='oversample').array)
    
    
    
class GaussianPSF(PSF):
    """
    Generate a PSF for an Gaussian function
    
    Needed keywords arguments:
    - fwhm: [arcsec] the FWHM of the PSF.
    
    Optional keywords
    - size: [int] the side length of the array in pixels
    - pix_res: [arcsec] the pixel scale used in the array, default is 0.004
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
        size = int(np.max((round(5 * self.fwhm / pix_res) * 2 + 1, size)))

        if size > 512: 
            size = 512
            print("FWHM [arcsec]:", fwhm, "- pixel res [arcsec]:", pix_res)
            print("Array size:", size,"x",size, "- PSF FoV:", size * pix_res)
            warnings.warn("PSF dimensions too large")
        
        super(GaussianPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Gaussian"
        self.info['description'] = "Gaussian PSF, FWHM = %.1f arcsec" \
                                    % (self.fwhm * 1E3)
                                                                  
        n = (self.fwhm / 2.35) / self.pix_res
        self.set_array(Gaussian2DKernel(n, x_size=self.size, y_size=self.size,
                                        mode='oversample').array)
        
        
        
class MoffatPSF(PSF):
    """
    Generate a PSF for a Moffat function. Alpha is generated from the FWHM and 
    Beta = 4.765 (from Trujillo et al. 2001)
    
    Needed keywords arguments:
    - fwhm: [arcsec] the FWHM of the PSF.
    
    Optional keywords
    - size: [int] the side length of the array in pixels
    - pix_res: [arcsec] the pixel scale used in the array, default is 0.004
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
        size = int(np.max((round(4 * self.fwhm / pix_res) * 2 + 1, size)))

        if size > 511: 
            size = 511
            print("FWHM [arcsec]:", fwhm, "- pixel res [arcsec]:", pix_res)
            print("Array size:", size,"x",size, "- PSF FoV:", size * pix_res)
            warnings.warn("PSF dimensions too large")
        
        super(MoffatPSF, self).__init__(size, pix_res)

        beta = 4.765 ### Trujillo et al. 2001
        alpha = self.fwhm/(2 * np.sqrt(2**(1/beta) - 1))
        self.info["Type"] = "Moffat"
        self.info['description'] = "Moffat PSF, FWHM = %.1f, alpha = %.1f"\
                                       % (self.fwhm * 1E3, alpha)

        if self.size > 100:
            mode = "linear_interp"
        else:
            mode = "oversample"
        self.set_array(Moffat2DKernel(alpha, beta, x_size=self.size, 
                                      y_size=self.size, mode=mode).array)
 
 

class CombinedPSF(PSF):
    """
    Generate a PSF from a collection of several PSFs. 
    
    Keywords:
    - psf_list: [list] A list of PSF objects
    
    Optional keywords:
    - size: [int] the side length in pixels of the array
    """ 

    def __init__(self, psf_list, **kwargs):
        """Generate a master psf through convolution of a list of psfs"""
        import copy
        
        if not hasattr(psf_list, "__len__") or len(psf_list) < 2:
            raise ValueError("psf_list requires more than 1 PSF object")

        pix_res_list = [psf.pix_res for psf in psf_list]
        if not all(res == pix_res_list[0] for res in pix_res_list):
            raise ValueError("Not all PSFs in have the same pixel resolution")

        pix_res = pix_res_list[0]
        
        if "size" in kwargs.keys():
            size = int(kwargs["size"] // 2) * 2 + 1
        else: 
            size_list = [psf.size for psf in psf_list]
            size = int(np.max(size_list) // 2) * 2 + 1
        
        ## Compensate for the shift in centre due to a DeltaPSF
        shifts = np.asarray([(0,0)] + [psf.position for psf in psf_list \
                                            if psf.info["Type"] == "Delta"])
        size += 2 * np.max(shifts)
        print(size)
                                             
        arr_tmp = np.zeros((size, size)) 
        arr_tmp[size // 2, size // 2] = 1
        
        for psf in psf_list:
            arr_tmp = convolve_fft(arr_tmp, psf.array)

                
        super(CombinedPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Combined"
        self.info['description'] = "Combined PSF from " + str(len(psf_list)) \
                                                                + "PSF objects"
        self.set_array(arr_tmp)

        
        
class UserPSF(PSF):
    """
    Import a PSF from a FITS file. 
    
    Keywords:
    - filename: [str] path to the FITS file to be read in
    
    Optional keywords
    - fits_ext: [int] the FITS extension number (default 0) for the data
    - pix_res: [arcsec] the pixel scale used in the array, default is 0.004
    """ 

    def __init__(self, filename, **kwargs):

        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]  
        else: 
            pix_res = 0.004
    
        if "fits_ext" in kwargs.keys():
            fits_ext = kwargs["fits_ext"]  
        else: 
            fits_ext = 0
            
        self.filename = filename
        self.fits_ext = fits_ext
        
        header = fits.getheader(self.filename, ext = self.fits_ext)
        data = fits.getdata(self.filename, ext = self.fits_ext)
        size = header["NAXIS1"]

        super(UserPSF, self).__init__(size, pix_res)
        self.info["Type"] = "User"
        self.info['description'] = "PSF from FITS file: " + self.filename
        
        self.set_array(data)

        if "size" in kwargs.keys():
            self.resize(kwargs["size"])


            
###############################################################################
#                       PSFCube and PSFCube subclasses                        #
###############################################################################
                 
# (Sub)Class    Needed              Optional
# PSF           size, pix_res
# DeltaPSF                          pix_res=0.004, size=f(position), position=(0,0)
# AiryPSF       fwhm                pix_res=0.004, size=f(fwhm)
# GaussianPSF   fwhm                pix_res=0.004, size=f(fwhm)
# MoffatPSF     fwhm                pix_res=0.004, size=f(fwhm)
# CombinedPSF   psf_list            size=f(psf_list)
# UserPSF       filename,           pix_res=0.004, size=f(filename), fits_ext=0

#    Keywords:
#    - type:
#    - lam_bin_centers
    
#    Bound keywords:
#    - fwhm : needed for AiryPSF, GaussianPSF, MoffatPSF
#    - psf_list: needed for CombinedPSF
#    - filename: needed for UserPSF
    
#    Optional keywords
#    - size:
#    - pix_res:
#    - position: optional in DeltaPSF
#    - fits_ext: optional in UserPSF


class PSFCube(object):
    """Class holding wavelength dependent point spread function

    Keywords:
    - lam_bin_centers: [µm] the centre of each wavelength slice
    """
    
    def __init__(self, lam_bin_centers):
            
        self.lam_bin_centers = lam_bin_centers
        self.psf_slices = [None] * len(lam_bin_centers)
        
        self.info = dict([])
        self.info['created'] = 'yes'
        self.info['description'] = "Point spread function (multiple layer)"
    
    def __repr__(self):
        return self.info['description']
        
    def __getitem__(self, i):
        return self.psf_slices[i]
        
    #def __array__(self):
     #   return self.psf_slices
        
    def __len__(self):
        return len(self.psf_slices)

    def resize(self, new_size):
        """Resize the list of PSFs. The target shape is (new_size, new_size).

        Keywords:
        - new_size: [int] the new size of the PSF array in pixels
        """
        for psf in psf_slices: 
            psf.resize(new_size)
    
    def resample(new_pix_res):
        """
        Resample the the list of PSF array onto a new grid and return the new
        grid - Not perfect, but conserves flux
        Example: new_PSF = old_PSF.resample(new_pix_res)
        
        Keywords:
        - new_pix_res: [arcsec] the pixel resolution of the returned array
        """
        return [psf.resample(new_pix_res) for psf in psf_slices]

    def export_to_FITS(self, filename, clobber=True, **header_info):
        """
        Export the PSFCube to a FITS file for later use
        
        Keywords:
        - filename
        - **header_info: A dict with extra header keys and values
        """
        ext_list = []

        for i in range(len(self)):
            psf = self.psf_slices[i]    
            if i == 0: 
                hdu = fits.PrimaryHDU(psf.array)
            else:
                hdu = fits.ImageHDU(psf.array)
            hdu.header["CDELT1"] = (psf.pix_res, "[arcsec] - Pixel resolution")
            hdu.header["CDELT2"] = (psf.pix_res, "[arcsec] - Pixel resolution")
            hdu.header["WAVECENT"] = (self.lam_bin_centers[i], "[micron] - Wavelength of slice")
            #hdu.header["PSF_TYPE"] = (self.info["Type"], "Type of PSF")
            ext_list += [hdu]
        
        hdu_list = fits.HDUList(ext_list)
        hdu_list.writeto(filename, clobber=clobber)
    

    def convolve(kernel_list):
        if len(self.psf_slices) != len(kernel_list):
            print("len(PSF_slices):",  len(self.psf_slices), 
                  "len(kernel_list):", len(kernel_list))
            raise ValueError("Number of kernels must equal number of PSFs")
        
        for psf, kernel in zip(self.psf_slices, kernel_list):
            psf.convolve(kernel)
          
    def __mul__(x):
        if not hasattr(x, "__len__"):
            y = [x] * len(self.psf_slices)
        else: 
            y = x
        
        if len(self.psf_slices) != len(y):
            print(len(self.psf_slices), len(y))
            raise ValueError("len(arguments) must equal len(PSFs)")
        
        for psf, y in zip(self.psf_slices, y):
                psf = psf * y

    def __add__(x):
        if not hasattr(x, "__len__"):
            y = [x] * len(self.psf_slices)
        else: 
            y = x
        
        if len(self.psf_slices) != len(y):
            print(len(self.psf_slices), len(y))
            raise ValueError("len(arguments) must equal len(PSFs)")
        
        for psf, y in zip(self.psf_slices, y):
                psf = psf + y
   
    def __sub__(x):
        if not hasattr(x, "__len__"):
            y = [x] * len(self.psf_slices)
        else: 
            y = x
        
        if len(self.psf_slices) != len(y):
            print(len(self.psf_slices), len(y))
            raise ValueError("len(arguments) must equal len(PSFs)")
        
        for psf, y in zip(self.psf_slices, y):
                psf = psf - y
    
    def __rmul__(self, x):
        self.__mul__(x)
                
    def __radd__(self, x):
        self.__add__(x)
    
    def __rsub__(self, x):
        self.__sub__(x)    


        
class DeltaPSFCube(PSFCube):
    """
    Generate a list of DeltaPSFs for wavelengths defined in lam_bin_centers
    
    Keywords:
    - lam_bin_centers: [µm] the centre of each wavelength slice
    
    Optional keywords:
    - positions: [(px,px)] either a tuple, or a list of tuples denoting the 
                 position of the delta function
    """

    def __init__(self, lam_bin_centers, positions=(0,0), **kwargs):
        super(DeltaPSFCube, self).__init__(lam_bin_centers)
        
        if not hasattr(positions[0], "__len__"):
            positions = [positions]*len(self), 
        
        for i in range(len(self)):
            self.psf_slices[i] = DeltaPSF(position=positions[i], **kwargs)
        
        self.info['description'] = "List of Delta function PSFs"
        self.info["Type"] = "DeltaCube"
    
class AiryPSFCube(PSFCube):
    """
    Generate a list of AiryPSFs for wavelengths defined in lam_bin_centers
    
    Keywords:
    - lam_bin_centers: [µm] a list with the centres of each wavelength slice
    
    Optional keywords:
    - fwhm: [arcsec] the equivalent FWHM of the PSF.
    - diameter: [m] diamter of primary mirror. Default is 39.3m.
    """
    
    def __init__(self, lam_bin_centers, fwhm=None, **kwargs):
        super(AiryPSFCube, self).__init__(lam_bin_centers)
        
        if "diameter" in kwargs.keys(): 
            self.diameter = kwargs["diameter"] 
        else:
            self.diameter = 39.3
            
        if fwhm is None:
            # lam in µm, diameter in m, 206265 is 1 rad in arcsec
            self.fwhm = [206265 * 1.22 * lam * 1E-6 / self.diameter \
                                                    for lam in lam_bin_centers]
        elif not hasattr(fwhm, "__len__"):
            self.fwhm = [fwhm] * len(self)
        
        self.psf_slices = [AiryPSF(fwhm = f, **kwargs) for f in self.fwhm]
    
        self.info['description'] = "List of Airy function PSFs"
        self.info["Type"] = "AiryCube"
    
    
class GaussianPSFCube(PSFCube):
    """
    Generate a list of GaussianPSFs for wavelengths defined in lam_bin_centers
    
    Keywords:
    - lam_bin_centers: [µm] a list with the centres of each wavelength slice 
    
    Optional keywords:
    - fwhm: [arcsec] the FWHM of the PSF.
   - diameter: [m] diamter of primary mirror. Default is 39.3m.
    """
    
    def __init__(self, lam_bin_centers, fwhm=None, **kwargs):
        super(GaussianPSFCube, self).__init__(lam_bin_centers)
        
        if "diameter" in kwargs.keys(): 
            self.diameter = kwargs["diameter"] 
        else:
            self.diameter = 39.3
            
        if fwhm is None:
            # lam in µm, diameter in m, 206265 is 1 rad in arcsec
            self.fwhm = [206265 * 1.22 * lam * 1E-6 / self.diameter \
                                                    for lam in lam_bin_centers]
        elif not hasattr(fwhm, "__len__"):
            self.fwhm = [fwhm] * len(self)
        
        self.psf_slices = [GaussianPSF(fwhm = f, **kwargs) for f in self.fwhm]
    
        self.info['description'] = "List of Gaussian function PSFs"
        self.info["Type"] = "GaussianCube"
    
class MoffatPSFCube(PSFCube):
    """
    Generate a list of MoffatPSFs for wavelengths defined in lam_bin_centers
    
    Keywords:
    - lam_bin_centers: [µm] a list with the centres of each wavelength slice 
    
    Optional keywords:
    - fwhm: [arcsec] the FWHM of the PSF.
    - diameter: [m] diamter of primary mirror. Default is 39.3m.
    """
    
    def __init__(self, lam_bin_centers, fwhm=None, **kwargs):
        super(MoffatPSFCube, self).__init__(lam_bin_centers)
        
        if "diameter" in kwargs.keys(): 
            self.diameter = kwargs["diameter"] 
        else:
            self.diameter = 39.3
            
        if fwhm is None:
            # lam in µm, diameter in m, 206265 is 1 rad in arcsec
            self.fwhm = [206265 * 1.22 * lam * 1E-6 / self.diameter \
                                                    for lam in lam_bin_centers]
        elif not hasattr(fwhm, "__len__"):
            self.fwhm = [fwhm] * len(self)
        
        self.psf_slices = [MoffatPSF(fwhm = f, **kwargs) for f in fwhm]
        
        self.info['description'] = "List of Moffat function PSFs"
        self.info["Type"] = "MoffatCube"
        
class CombinedPSFCube(PSFCube):
    """
    Generate a list of CombinedPSFs from the list of PSFCubes in psfcube_list
    
    Keywords:
    - lam_bin_centers: [µm] a list with the centres of each wavelength slice 
    
    Optional keywords:
    - fwhm: [arcsec] the FWHM of the PSF.
    - diameter: [m] diamter of primary mirror. Default is 39.3m.
    """
    
    def __init__(self, psfcube_list, **kwargs):
        
        if not (type(psfcube_list) == list and len(psfcube_list) >= 2):
            raise ValueError("psfcube_list only takes a list of PSFCube objects")
            
        ## Check that the wavelengths are equal
        lam_list = [cube.lam_bin_centers for cube in psfcube_list]
        if not all([all(lam == lam_list[0]) for lam in lam_list]):
            raise ValueError("Wavelength arrays of psf cubes are not equal")
        lam_bin_centers = lam_list[0]

        super(CombinedPSFCube, self).__init__(lam_bin_centers)

        self.info['description'] = "Master psf cube from list"
        for i in range(len(psfcube_list)):
            self.info['PSF%02d' % (i+1)] = psfcube_list[i].info['description']

            
        for i in range(len(self)):
            self.psf_slices[i] = CombinedPSF([psf[i] for psf in psfcube_list], **kwargs)
            
        
    
class UserPSFCube(PSFCube):
    """
    """
    

    
class ADC_PSFCube(PSFCube):
    """
    """
    
 
    
    
    
    
    
    
    
    
    
    
    
    
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