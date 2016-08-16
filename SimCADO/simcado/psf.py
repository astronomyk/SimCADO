"""psf.py"""
###############################################################################
# PSFs and PSFCubes
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
#   - UserPSFCube        (filename, ext_no=0)
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
# A psf instance will have these additional arguments:
#   - cube ... a (l,x,y) 3D array to hold the PSF cube
#
# As far as input goes, psf should be able to accept a dictionary with the
# keywords necessary to build the cube.
#
# NOTES:
# All wavelength values are given in [um]
# All pixel dimensions are given in [arcsec]
# All angles are given in [deg]
#
#
# CLASSES:
#  PSF(object)
#  psf(object)
#
#
# SUBCLASSES:
#  MoffatPSF(PSF)
#  GaussianPSF(PSF)
#  AiryPSF(PSF)
#  DeltaPSF(PSF)
#  LinePSF(PSF)
#  UserPSFCube(PSF)
#
#  Deltapsf(psf)
#  Airypsf(psf)
#  Gaussianpsf(psf)
#  Moffatpsf(psf)
#  CombinedPSFCube(psf)
#  UserPSFCube(psf)
#  ADC_psf(psf)
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
# Both PSF and psf can be created from a single model or through
# convolution of a list of PSF components

import warnings
from copy import deepcopy

import numpy as np
import scipy.ndimage.interpolation as spi

from astropy.io import fits
#from astropy import units as u   ## unused (OC)
from astropy.convolution import (Gaussian2DKernel,
                                 Moffat2DKernel)
from astropy.convolution import convolve_fft
from astropy.convolution import Kernel2D
from astropy.modeling.core import Fittable2DModel
from astropy.modeling.parameters import Parameter

try:
    import simcado.utils as utils
except ImportError:
    import utils



## TODO
# - Add a ellipticity to the GaussianPSF
# - Make sure MoffatPSF works



## These classes and functions are exported to the package
__all__ = ["PSF", "PSFCube",
           "MoffatPSF", "MoffatPSFCube",
           "AiryPSF", "AiryPSFCube",
           "GaussianPSF", "GaussianPSFCube",
           "DeltaPSF", "DeltaPSFCube",
           "CombinedPSF", "CombinedPSFCube",
           "UserPSF", "UserPSFCube"]


###############################################################################
#                            PSF and PSF subclasses                           #
###############################################################################

# (Sub)Class    Needed              Optional
# PSF           size, pix_res
# DeltaPSF                          pix_res=0.004, size=f(position), position=(0,0)
# AiryPSF       fwhm                pix_res=0.004, size=f(fwhm)
# GaussianPSF   fwhm                pix_res=0.004, size=f(fwhm)
# MoffatPSF     fwhm                pix_res=0.004, size=f(fwhm)
# CombinedPSFCube   psf_list            size=f(psf_list)
# UserPSFCube       filename,           pix_res=0.004, size=f(filename), fits_ext=0



class PSF(object):
    """Point spread function (single layer) base class


    Parameters
    ==========
    - size : int
        [pixel] the side length of the array
    - pix_res : float
        [arcsec] the pixel scale of the array
    """

    def __init__(self, size, pix_res):
        self.size = size
        self.shape = [size, size]
        self.pix_res = pix_res
        self.array = np.zeros((self.size, self.size))

        self.info = dict([])
        self.info['description'] = "Point spread function (single layer)"

    def __str__(self):
        return self.info['description']

    def set_array(self, array, threshold=1e-15):
        """
        Renormalise the array and make sure there aren't any negative values
        that will screw up the flux later on

        Parameters
        ==========
        - array: [ndarray] the array representing the PSF
        - threshold: by default set to 1E-15
        """
        self.array = array.astype(np.float32)
        self.array[self.array <= 0] = threshold
        self.array = self.array / np.sum(self.array)
        self.size = self.array.shape[0]
        self.shape = self.array.shape

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
        new_arr *= np.sum(self.array) / np.sum(new_arr)

        # catch any even-numbered arrays
        if new_arr.shape[0] % 2 == 0:
            tmp_arr = np.zeros((new_arr.shape[0]+1, new_arr.shape[1]+1))
            tmp_arr[:new_arr.shape[0], :new_arr.shape[1]] = new_arr
            new_arr = spi.shift(tmp_arr, 0.5, order=5)

        ############################################################
        # Not happy with the way the returned type is not the same #
        # as the original type. The new object is a plain PSF      #
        ############################################################
        new_psf = PSF(size=new_arr.shape[0], pix_res=new_pix_res)
        new_psf.set_array(new_arr)
        return new_psf

    def convolve(self, kernel):
        """
        Convolve the PSF with another kernel. The PSF keeps its shape

        Keywords:
        - kernel: [PSF/ndarray] either a numpy.ndarray or a PSF (sub)class
        """
        psf_new = deepcopy(self)

        if issubclass(type(kernel), PSF):
            psf_new.set_array(convolve_fft(self.array, kernel.array))
        else:
            psf_new.set_array(convolve_fft(self.array, kernel))
        psf_new.info["Type"] = "Combined"

        return psf_new

    def __array__(self):
        return self.array

    def __mul__(self, x):
        psf_new = deepcopy(self)
        return psf_new.array * x

    def __add__(self, x):
        psf_new = deepcopy(self)
        return psf_new.array + x

    def __sub__(self, x):
        psf_new = deepcopy(self)
        return psf_new.array - x

    def __rmul__(self, x):
        return self.__mul__(x)

    def __radd__(self, x):
        return self.__add__(x)

    def __rsub__(self, x):
        psf_new = deepcopy(self)
        return x - psf_new.array

    def __imul__(self, x):
        return self.__mul__(x)

    def __iadd__(self, x):
        return self.__add__(x)

    def __isub__(self, x):
        return self.__sub__(x)


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
            self.position = (0, 0)

        if "size" in kwargs.keys():
            size = round(kwargs["size"] / 2) * 2 + 5
        else:
            size = int(np.max(np.abs(self.position))) * 2 + 5

        if not np.max(self.position) < size:
            raise ValueError("positions are outside array borders:")

        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]
        else:
            pix_res = 0.004

        super(DeltaPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Delta"
        self.info["description"] = "Delta PSF, centred at (%.1f, %.1f)" \
                                    % self.position
        self.info["x_shift"] = self.position[0]
        self.info["y_shift"] = self.position[1]

        self.x = self.size // 2 + self.position[0]
        self.y = self.size // 2 + self.position[1]

        x2 = self.x - int(self.x)
        x1 = 1. - x2
        y2 = self.y - int(self.y)
        y1 = 1. - y2

        arr = np.zeros((self.size, self.size))
        arr[int(self.y) : int(self.y) + 2, int(self.x) : int(self.x) + 2] = \
            np.array([[x1 * y1, x2 * y1], [x1 * y2, x2 * y2]])
        self.set_array(arr)



class AiryPSF(PSF):
    """
    Generate a PSF for an Airy function with an equivalent FWHM

    Needed keywords arguments:
    - fwhm: [arcsec] the equivalent FWHM of the Airy disk core.

    Optional keywords
    - size [int] the side length of the array in pixels
    - pix_res [arcsec] the pixel scale used in the array, default is 0.004
    - obscuration [0..1]: radius of inner obscuration as fraction of aperture
      radius
    - mode: ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)
    """

    def __init__(self, fwhm, obscuration=0., size=255, pix_res=0.004,
                 **kwargs):

        # Ensure that PSF size is at least 8 times fwhm
        size = int(np.max((round(8 * fwhm / pix_res) * 2 + 1, size)))

        # if size > 511:
            # size = 511
            # print("FWHM [arcsec]:", fwhm, "- pixel res [arcsec]:", pix_res)
            # print("Array size:", size, "x", size, "- PSF FoV:", size * pix_res)
            # warnings.warn("PSF dimensions too large - cropped to 512x512")

        # Check for 'mode' keyword argument
        if "mode" in kwargs.keys():
            mode = kwargs["mode"]
        else:
            if size > 100:
                mode = "linear_interp"
            else:
                mode = 'oversample'

        # Ensure that obscuration is between 0 and 1
        if obscuration < 0 or obscuration > 1:
            print("Warning: Obscuration must be between 0 and 1. Using 0.")
            obscuration = 0.

        super(AiryPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Airy"
        self.info['description'] \
            = "Airy PSF, FWHM = %.1f mas, obscuration = %.2f" \
            % (fwhm * 1E3, obscuration)
        self.info["fwhm"] = fwhm * 1E3           # milliarcseconds
        self.info["obscuration"] = obscuration
        self.info["pixel scale"] = pix_res
        self.info["mode"] = mode

        self.fwhm = fwhm
        self.pix_res = pix_res

        # These are first zero and FWHM of J1(x)/x.
        airy_zero = 3.8317059860286755
        airy_fwhm = 2 * 1.616339948310703

        # rescale radial array so that the PSF has the required FWHM
        # (for eps = 0)
        first_zero = fwhm / airy_fwhm * airy_zero / pix_res
        psf_arr = AiryDiskDiff2DKernel(first_zero,
                                       obscuration=obscuration,
                                       x_size=size,
                                       y_size=size,
                                       mode=mode).array
        self.set_array(psf_arr)



class GaussianPSF(PSF):
    """
    Generate a PSF for an Gaussian function

    Needed keywords arguments:
    - fwhm: [arcsec] the FWHM of the PSF.

    Optional keywords
    - size: [int] the side length of the array in pixels
    - pix_res: [arcsec] the pixel scale used in the array, default is 0.004
    - mode: ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)
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

        if "undersized" in kwargs.keys() and kwargs["undersized"]:
            pass
        else:
            size = int(np.max((round(5 * self.fwhm / pix_res) * 2 + 1, size)))

        # Check for 'mode' keyword argument
        if "mode" in kwargs.keys():
            mode = kwargs["mode"]
        else:
            if size > 100:
                mode = "linear_interp"
            else:
                mode = 'oversample'

        super(GaussianPSF, self).__init__(size, pix_res)
        self.info["Type"] = "Gaussian"
        self.info['description'] = "Gaussian PSF, FWHM = %.1f mas" \
                                    % (self.fwhm * 1E3)
        self.info["fwhm"] = self.fwhm * 1E3

        n = (self.fwhm / 2.35) / self.pix_res
        k = Gaussian2DKernel(n, x_size=self.size, y_size=self.size, mode=mode).array

        self.set_array(k)



class MoffatPSF(PSF):
    """
    Generate a PSF for a Moffat function. Alpha is generated from the FWHM and
    Beta = 4.765 (from Trujillo et al. 2001)

    Needed keywords arguments:
    - fwhm: [arcsec] the FWHM of the PSF.

    Optional keywords
    - size: [int] the side length of the array in pixels
    - pix_res: [arcsec] the pixel scale used in the array, default is 0.004
    - mode: ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)
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

        # Check for 'mode' keyword argument
        if "mode" in kwargs.keys():
            mode = kwargs["mode"]
        else:
            if size > 100:
                mode = "linear_interp"
            else:
                mode = 'oversample'

        super(MoffatPSF, self).__init__(size, pix_res)

        beta = 4.765 ### Trujillo et al. 2001
        alpha = self.fwhm/(2 * np.sqrt(2**(1/beta) - 1))
        self.info["Type"] = "Moffat"
        self.info['description'] = "Moffat PSF, FWHM = %.1f, alpha = %.1f"\
                                       % (self.fwhm * 1E3, alpha)
        self.info["fwhm"] = self.fwhm * 1E3

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

        if not hasattr(psf_list, "__len__") or len(psf_list) < 2:
            raise ValueError("psf_list requires more than 1 PSF object")

        pix_res_list = [psf.pix_res for psf in psf_list]
        if not all(res == pix_res_list[0] for res in pix_res_list):
            raise ValueError("Not all PSFs have the same pixel resolution")

        pix_res = pix_res_list[0]

        if "size" in kwargs.keys():
            size = int(kwargs["size"] // 2) * 2 + 1
        else:
            size_list = [psf.size for psf in psf_list]
            size = int(np.max(size_list) // 2) * 2 + 1

        ## Compensate for the shift in centre due to a DeltaPSF
        shifts = np.asarray([(0, 0)] + [psf.position for psf in psf_list \
                                        if psf.info["Type"] == "Delta"])
        size += 2 * np.max(shifts)

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

        if "fits_ext" in kwargs.keys():
            fits_ext = kwargs["fits_ext"]
        else:
            fits_ext = 0

        self.filename = filename
        self.fits_ext = fits_ext

        header = fits.getheader(self.filename, ext=self.fits_ext)
        data = fits.getdata(self.filename, ext=self.fits_ext)
        size = header["NAXIS1"]

        if "pix_res" in kwargs.keys():
            pix_res = kwargs["pix_res"]
        elif "CDELT1" in header.keys():
            pix_res = header["CDELT1"]
        else:
            pix_res = 0.004

        super(UserPSF, self).__init__(size, pix_res)
        self.info["Type"] = "User"
        self.info['description'] = "PSF from FITS file: " + self.filename

        self.set_array(data)

        if "size" in kwargs.keys():
            self.resize(kwargs["size"])



###############################################################################
#                       psf and psf subclasses                        #
###############################################################################

# (Sub)Class    Needed              Optional
# PSF           size, pix_res
# DeltaPSF                          pix_res=0.004, size=f(position), position=(0,0)
# AiryPSF       fwhm                pix_res=0.004, size=f(fwhm)
# GaussianPSF   fwhm                pix_res=0.004, size=f(fwhm)
# MoffatPSF     fwhm                pix_res=0.004, size=f(fwhm)
# CombinedPSFCube   psf_list            size=f(psf_list)
# UserPSFCube       filename,           pix_res=0.004, size=f(filename), fits_ext=0

#    Keywords:
#    - type:
#    - lam_bin_centers

#    Bound keywords:
#    - fwhm : needed for AiryPSF, GaussianPSF, MoffatPSF
#    - psf_list: needed for CombinedPSFCube
#    - filename: needed for UserPSFCube

#    Optional keywords
#    - size:
#    - pix_res:
#    - position: optional in DeltaPSF
#    - fits_ext: optional in UserPSFCube


class PSFCube(object):
    """Class holding wavelength dependent point spread function.
    Special functions:
    - len(self) return the number of layers in the psf
    - self[i] returns the PSF object for layer i. If the __array__ function
      is called, self[i] will return the array associated with the instance
      e.g plt.imshow(self[i]) will plot PSF.array from self.psf_slices[i]
    - Maths operators *,+,- act equally on all PSF.arrays in self.psf_slices

    Keywords:
    - lam_bin_centers: [um] the centre of each wavelength slice
    """

    def __init__(self, lam_bin_centers):

        self.lam_bin_centers = lam_bin_centers
        self.psf_slices = [None] * len(lam_bin_centers)

        self.info = dict([])
        self.info['created'] = 'yes'
        self.info['description'] = "Point spread function (multiple layer)"


    def resize(self, new_size):
        """Resize the list of PSFs. The target shape is (new_size, new_size).

        Keywords:
        - new_size: [int] the new size of the PSF array in pixels
        """
        for psf in self.psf_slices:
            psf.resize(new_size)

    def resample(self, new_pix_res):
        """
        Resample the the list of PSF array onto a new grid and return the new
        grid - Not perfect, but conserves flux
        Example: new_PSF = old_PSF.resample(new_pix_res)

        Keywords:
        - new_pix_res: [arcsec] the pixel resolution of the returned array
        """
        ## TODO: Check whether this makes sense
        self.psf_slices = [psf.resample(new_pix_res) for psf in self.psf_slices]

    def export_to_fits(self, filename, clobber=True, **header_info):
        """
        Export the psf to a FITS file for later use

        Keywords:
        - filename
        - **header_info: A dict with extra header keys and values
        """
        ## TODO: use header_info or drop it (OC)

        ext_list = fits.HDUList()

        #for i in range(len(self)):
        #    psf = self.psf_slices[i]
        for i, psf in enumerate(self.psf_slices):
            #if i == 0:
            #    hdu = fits.PrimaryHDU(psf.array)
            #else:
            ## PrimaryHDUs are typically empty. Not necessary (OC)
            hdu = fits.ImageHDU(psf.array)
            hdu.header["CDELT1"] = (psf.pix_res, "[arcsec] - Pixel resolution")
            hdu.header["CDELT2"] = (psf.pix_res, "[arcsec] - Pixel resolution")
            hdu.header["WAVE0"] = (self.lam_bin_centers[i],
                                   "[micron] - Wavelength of slice")
            hdu.header["NSLICES"] = (len(self), "Number of wavelength slices")

            for k in psf.info.keys():
                hdu.header[k[:8].upper()] = (self[i].info[k], k)

            ext_list.append(hdu)

        ext_list.writeto(filename, clobber=clobber)


    def convolve(self, kernel_list):
        """Convolve a list of PSFs with a list of kernels"""
        if len(self.psf_slices) != len(kernel_list):
            print("len(PSF_slices):", len(self.psf_slices),
                  "len(kernel_list):", len(kernel_list))
            raise ValueError("Number of kernels must equal number of PSFs")

        for i in np.arange(len(self.psf_slices)):
            tmp = self.psf_slices[i].convolve(kernel_list[i])
            self.psf_slices[i].set_array(tmp.array)
        self.info["Type"] = "Complex"


    def nearest(self, lam):
        """Returns the PSF closest to the desired wavelength, lam [um] """
        i = utils.nearest(self.lam_bin_centers, lam)
        return self.psf_slices[i]


    def __str__(self):
        return self.info['description']

    def __getitem__(self, i):
        if len(self.psf_slices) > 1:
            return self.psf_slices[i]
        else:
            return self.psf_slices[0]

    def __len__(self):
        return len(self.psf_slices)


    def __mul__(self, x):
        newpsf = deepcopy(self)

        if not hasattr(x, "__len__"):
            y = [x] * len(self.psf_slices)
        else:
            y = x

        if len(self.psf_slices) != len(y):
            print(len(self.psf_slices), len(y))
            raise ValueError("len(arguments) must equal len(PSFs)")

        for i in np.arange(len(self.psf_slices)):
            newpsf[i].set_array(self.psf_slices[i] * y[i])
        return newpsf



    def __add__(self, x):
        newpsf = deepcopy(self)

        if not hasattr(x, "__len__"):
            y = [x] * len(self.psf_slices)
        else:
            y = x

        if len(self.psf_slices) != len(y):
            print(len(self.psf_slices), len(y))
            raise ValueError("len(arguments) must equal len(PSFs)")

        for i in np.arange(len(self.psf_slices)):
            newpsf[i].set_array(self.psf_slices[i] + y[i])
        return newpsf


    def __sub__(self, x):
        newpsf = deepcopy(self)

        if not hasattr(x, "__len__"):
            y = [x] * len(self.psf_slices)
        else:
            y = x

        if len(self.psf_slices) != len(y):
            print(len(self.psf_slices), len(y))
            raise ValueError("len(arguments) must equal len(PSFs)")

        for i in np.arange(len(self.psf_slices)):
            newpsf[i].set_array(self.psf_slices[i] - y[i])
        return newpsf



    def __rmul__(self, x):
        self.__mul__(x)

    def __radd__(self, x):
        self.__add__(x)

    def __rsub__(self, x):
        self.__sub__(x)



class DeltaPSFCube(PSFCube):
    """
    Generate a list of DeltaPSFs for wavelengths defined in lam_bin_centers

    Parameters:
    ===========
    - lam_bin_centers: [um] the centre of each wavelength slice

    Optional Parameters:
    ====================
    - positions: [(x,y)] either a tuple, or a list of tuples denoting the
                 position of the delta function
    - pix_res: [arcsec], pixel scale of the PSF. Default is 0.004 arcsec
    - size: [int], side length of the PSF array
    """

    def __init__(self, lam_bin_centers, positions=(0, 0), **kwargs):
        super(DeltaPSFCube, self).__init__(lam_bin_centers)

        if not hasattr(positions[0], "__len__"):
            positions = [positions]*len(self)

        for i in range(len(self)):
            self.psf_slices[i] = DeltaPSF(position=positions[i], **kwargs)

        self.info['description'] = "List of Delta function PSFs"
        self.info["Type"] = "DeltaCube"


class AiryPSFCube(PSFCube):
    """
    Generate a list of AiryPSFs for wavelengths defined in lam_bin_centers

    Keywords:
    - lam_bin_centers: [um] a list with the centres of each wavelength slice

    Optional keywords:
    - fwhm: [arcsec] the equivalent FWHM of the PSF.
    - diameter: [m] diamter of primary mirror. Default is 39.3m.
    - mode: ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)
    """

    def __init__(self, lam_bin_centers, fwhm=None, **kwargs):
        super(AiryPSFCube, self).__init__(lam_bin_centers)

        if "diameter" in kwargs.keys():
            self.diameter = kwargs["diameter"]
        else:
            self.diameter = 39.3

        if "obscuration" in kwargs.keys():
            self.obscuration = kwargs["obscuration"]
        else:
            self.obscuration = 11.1/39.3

        if fwhm is None:
            # lam in um, diameter in m, 206265 is 1 rad in arcsec
            self.fwhm = [206265 * 1.22 * lam * 1E-6 / self.diameter \
                                                    for lam in lam_bin_centers]
        elif not hasattr(fwhm, "__len__"):
            self.fwhm = [fwhm] * len(self)
        else:
            self.fwhm = fwhm

        self.psf_slices = [AiryPSF(fwhm=f, **kwargs) for f in self.fwhm]
        self.size = [psf.size for psf in self.psf_slices]

        self.info['description'] = "List of Airy function PSFs"
        self.info["Type"] = "AiryCube"


class GaussianPSFCube(PSFCube):
    """
    Generate a list of GaussianPSFs for wavelengths defined in lam_bin_centers

    Keywords:
    - lam_bin_centers: [um] a list with the centres of each wavelength slice

    Optional keywords:
    - fwhm: [arcsec] the FWHM of the PSF.
    - diameter: [m] diamter of primary mirror. Default is 39.3m.
    - mode: ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)
    """

    def __init__(self, lam_bin_centers, fwhm=None, **kwargs):
        super(GaussianPSFCube, self).__init__(lam_bin_centers)

        if "diameter" in kwargs.keys():
            self.diameter = kwargs["diameter"]
        else:
            self.diameter = 39.3

        if fwhm is None:
            # lam in um, diameter in m, 206265 is 1 rad in arcsec
            self.fwhm = [206265 * 1.22 * lam * 1E-6 / self.diameter \
                                                    for lam in lam_bin_centers]
        elif not hasattr(fwhm, "__len__"):
            self.fwhm = [fwhm] * len(self)

        self.psf_slices = [GaussianPSF(fwhm=f, **kwargs) for f in self.fwhm]
        self.size = [psf.size for psf in self.psf_slices]

        self.info['description'] = "List of Gaussian function PSFs"
        self.info["Type"] = "GaussianCube"

class MoffatPSFCube(PSFCube):
    """
    Generate a list of MoffatPSFs for wavelengths defined in lam_bin_centers

    Keywords:
    - lam_bin_centers: [um] a list with the centres of each wavelength slice

    Optional keywords:
    - fwhm: [arcsec] the FWHM of the PSF.
    - diameter: [m] diamter of primary mirror. Default is 39.3m.
    - mode: ['oversample','linear_interp'] see Kernel2D (scipy.convolution.core)
    """

    def __init__(self, lam_bin_centers, fwhm=None, **kwargs):
        super(MoffatPSFCube, self).__init__(lam_bin_centers)

        if "diameter" in kwargs.keys():
            self.diameter = kwargs["diameter"]
        else:
            self.diameter = 39.3

        if fwhm is None:
            # lam in um, diameter in m, 206265 is 1 rad in arcsec
            rad2arcsec = 3600 * 180. / np.pi
            self.fwhm = rad2arcsec * 1.22 * lam_bin_centers * 1E-6 / self.diameter
        elif not hasattr(fwhm, "__len__"):
            self.fwhm = [fwhm] * len(self)

        self.psf_slices = [MoffatPSF(fwhm=f, **kwargs) for f in fwhm]
        self.size = [psf.size for psf in self.psf_slices]

        self.info['description'] = "List of Moffat function PSFs"
        self.info["Type"] = "MoffatCube"

class CombinedPSFCube(PSFCube):
    """
    Generate a list of CombinedPSFCubes from the list of psfs in psf_list

    Keywords:
    - lam_bin_centers: [um] a list with the centres of each wavelength slice

    Optional keywords:
    - fwhm: [arcsec] the FWHM of the PSF.
    - diameter: [m] diamter of primary mirror. Default is 39.3m.
    """

    def __init__(self, psf_list, **kwargs):

        if not (isinstance(psf_list, list) and len(psf_list) >= 2):
            raise ValueError("psf_list only takes a list of psf objects")

        ## Check that the wavelengths are equal
        lam_list = [cube.lam_bin_centers for cube in psf_list]
        if not all([all(lam == lam_list[0]) for lam in lam_list]):
            raise ValueError("Wavelength arrays of psf cubes are not equal")
        lam_bin_centers = lam_list[0]

        super(CombinedPSFCube, self).__init__(lam_bin_centers)

        self.info['description'] = "Master psf cube from list"
        self.info["Type"] = "CombinedCube"

        for i, psfi in enumerate(psf_list):
            self.info['PSF%02d' % (i+1)] = psfi.info['description']

        for i in range(len(self)):
            self.psf_slices[i] = CombinedPSFCube([psf[i] for psf in psf_list],
                                                 **kwargs)
        self.size = [psf.size for psf in self.psf_slices]


class UserPSFCube(PSFCube):
    """
    Read in a psf previously saved as a FITS file
    Keywords needed for a psf to be read in:
    NSLICES, WAVECENT, NAXIS1, CDELT1, PSF_TYPE, DESCRIPT

    N.B. A separate function will exist to convert foreign PSF FITS files into
    psf readable FITS files

    Keywords:
    - filename: the path to the FITS file holding the cube
    """

    def __init__(self, filename, lam_bin_centers):
        if not hasattr(lam_bin_centers, "__len__"):
            lam_bin_centers = [lam_bin_centers]

        n_slices = len(fits.open(filename).info(output=False))
        psf_slices = []

        # pull out the wavelengths in the PSF FITS files
        psf_lam_cen = []
        hdulist = fits.open(filename)
        for i in range(n_slices):
            if "WAVE0" in hdulist[i].header.keys():
                psf_lam_cen += [hdulist[i].header["WAVE0"]]
            elif "WAVELENG" in hdulist[i].header.keys():
                psf_lam_cen += [hdulist[i].header["WAVELENG"]]
            else:
                raise ValueError("""Could not determine wavelength of PSF in
                                 extension """ + str(i) + """. FITS file
                                 needs either WAVE0 or WAVELENG header
                                 keywords. \n Use SimCADO.utils.add_keyword()
                                 to add WAVELENG to the FITS header""")
            # If the wavelength is not in the 0.1-2.5 range, it must be in [m]
            if psf_lam_cen[i] < 0.1:
                psf_lam_cen[i] *= 1E6
        hdulist.close()

        # find the closest PSFs in the file to what is needed for the psf
        i_slices = utils.nearest(psf_lam_cen, lam_bin_centers)

        # import only the relevant PSFs
        for i in i_slices:

            hdr = fits.getheader(filename, ext=i)
            self.header = hdr

            pix_res = hdr["CDELT1"]
            if pix_res > 1:
                warnings.warn("CDELT > 1. Assuming the scale to be [mas]")
                pix_res *= 1E-3

            psf = PSF(size=hdr["NAXIS1"], pix_res=pix_res)
            psf.set_array(fits.getdata(filename, ext=i))

            if "PSF_TYPE" in hdr.keys():
                psf.info["Type"] = hdr["PSF_TYPE"]
            else:
                psf.info["Type"] = "Unknown"

            if "DESCRIPT" in hdr.keys():
                psf.info["description"] = hdr["DESCRIPT"]
            else:
                psf.info["description"] = "Unknown"

            psf_slices += [psf]

        lam_bin_centers = np.array(lam_bin_centers)

        super(UserPSFCube, self).__init__(lam_bin_centers)
        self.psf_slices = psf_slices
        self.size = [psf.size for psf in self.psf_slices]

        self.info['description'] = "User PSF cube input from " + filename
        self.info["Type"] = psf_slices[0].info["Type"]+"Cube"



class ADC_PSFCube(DeltaPSFCube):
    """
    Generates a Deltapsf with the shifts required to mimic the ADC at a
    certain efficiency

    Keywords:
    - lam_bin_centers: [um] a list with the centres of each wavelength slice

    Optional Keywords:
    - pix_res: [arcsec] the pixel scale used in the array, default is 0.004
    - OBS_PARALLACTIC_ANGLE: [deg] the orientation of the input cube relative to
      the zenith
    - INST_ADC_EFFICIENCY: [%] efficiency of the ADC
    - SCOPE_LATITUDE: [deg] latitude of the telescope site
    - SCOPE_ALTITUDE: [m] hight above sea level of the telescope site
    - ATMO_REL_HUMIDITY: [%] relative humidity in percent
    - OBS_ZENITH_DIST: [deg] zenith distance of the object
    - ATMO_TEMPERATURE: [deg C] air temperature of the observing site in Celsius
    - ATMO_PRESSURE: [mbar] air pressure of the observing site in millibar

    The default values for the above mentioned keywords are:
    0.004 arcsec, 0 deg, 100%, -24.5 deg, 3064m, 60%, 60 deg, 0 C, 750mbar
    """

    def __init__(self, lam_bin_centers, **kwargs):
        params = {"pix_res"           :0.004,
                  "PARALLACTIC_ANGLE"  :0,
                  "INST_ADC_EFFICIENCY":100,
                  "SCOPE_LATITUDE"     :-24.5,
                  "SCOPE_ALTITUDE"     :3064,
                  "ATMO_REL_HUMIDITY"  :60,
                  "OBS_ZENITH_DIST"    :60,
                  "ATMO_TEMPERATURE"   :0,
                  "ATMO_PRESSURE"      :750}

        params.update(**kwargs)
        pix_res = params["pix_res"]
        para_angle = params["OBS_PARALLACTIC_ANGLE"]
        effectiveness = params["INST_ADC_PERFORMANCE"] / 100.

        ## get the angle shift for each slice
        angle_shift = [utils.atmospheric_refraction(lam,
                                                    params["OBS_ZENITH_DIST"],
                                                    params["ATMO_TEMPERATURE"],
                                                    params["ATMO_REL_HUMIDITY"],
                                                    params["ATMO_PRESSURE"],
                                                    params["SCOPE_LATITUDE"],
                                                    params["SCOPE_ALTITUDE"])
                       for lam in lam_bin_centers]

        ## convert angle shift into number of pixels
        ## pixel shifts are defined with respect to last slice
        pixel_shift = (angle_shift - angle_shift[-1]) / pix_res
        if np.max(np.abs(pixel_shift)) > 1000:
            raise ValueError("Pixel shifts too great (>1000), check units")

        ## Rotate by the paralytic angle
        x = -pixel_shift * np.sin(np.deg2rad(para_angle)) * (1. - effectiveness)
        y = -pixel_shift * np.cos(np.deg2rad(para_angle)) * (1. - effectiveness)
        positions = [(xi, yi) for xi, yi in zip(x, y)]

        super(ADC_PSFCube, self).__init__(lam_bin_centers,
                                          positions=positions,
                                          pix_res=pix_res)
        self.info["Type"] = "ADC_psf"
        self.info['description'] = "ADC PSF cube for ADC effectiveness:" + \
                                    str(params["INST_ADC_EFFICIENCY"]) + \
                                    ", z0:" + str(params["OBS_ZENITH_DIST"])




## The following two classes implement a kernel for the PSF of a centrally
## obscured circular aperture. The classes are modelled after the kernels
## in astropy.convolution.kernel and the models in astropy.modeling.models,
## both from astropy version 1.1.1.
class AiryDiskDiff2DKernel(Kernel2D):
    """
    2D kernel for PSF for annular aperture

    This kernel models the diffraction pattern of a circular aperture with
    a central circular obscuration. This kernel is normalized to a peak
    value of 1.

    Parameters
    ----------
    radius : float
        The radius of the unobscured Airy disk kernel (radius of the first
        zero). Compute this from the outer aperture radius.
    obscuration : float
        Fraction of the aperture that is obscured (inner radius / outer radius)
        Default obscuration = 0.
    x_size : odd int, optional
        Size in x direction of the kernel array. Default = 8 * radius.
    y_size : odd int, optional
        Size in y direction of the kernel array. Default = 8 * radius.
    mode : str, optional
        One of the following discretization modes:
            * 'center' (default)
                Discretize model by taking the value
                at the center of the bin.
            * 'linear_interp'
                Discretize model by performing a bilinear interpolation
                between the values at the corners of the bin.
            * 'oversample'
                Discretize model by taking the average
                on an oversampled grid.
            * 'integrate'
                Discretize model by integrating the
                model over the bin.
    factor : number, optional
        Factor of oversampling. Default factor = 10.

    See Also
    --------
    Gaussian2DKernel, Box2DKernel, Tophat2DKernel, MexicanHat2DKernel,
    Ring2DKernel, TrapezoidDisk2DKernel, AiryDisk2DKernel, Moffat2DKernel
    (in astropy.kernels)

    Examples
    --------
    Kernel response:

     .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from SimCADO.psf import AiryDiskDiff2DKernel
        airydiskdiff_2D_kernel = AiryDiskDiff2DKernel(10)
        plt.imshow(airydiskdiff_2D_kernel, interpolation='none', origin='lower')
        plt.xlabel('x [pixels]')
        plt.ylabel('y [pixels]')
        plt.colorbar()
        plt.show()
    """
    _is_bool = False


    def __init__(self, radius, obscuration, **kwargs):
        from astropy.convolution.kernels import _round_up_to_odd_integer

        self._model = AiryDiskDiff2D(1, 0, 0, radius, obscuration)
        self._default_size = _round_up_to_odd_integer(8 * radius)
        super(AiryDiskDiff2DKernel, self).__init__(**kwargs)
        self.normalize()
        self._truncation = None

class AiryDiskDiff2D(Fittable2DModel):
    """
    Two dimensional Airy disk model with central obscuration.

    Parameters
    ----------
    amplitude : float
        Amplitude of the function.
    x_0 : float
        x position of the maximum of the function
    y_0 : float
        y position of the maximum of the function
    radius : float
        The radius of the unobscured Airy disk (radius of the first zero).
    eps : float
        The ratio of the inner to the outer radius of an annular aperture.

    See Also
    --------
    AiryDisk2D

    Notes
    -----
    Model formula:

        .. math:: f(r) = A \\left[\\frac{2 J_1(\\frac{\\pi r}{R/R_z})}\\right]^2

    Where :math: `J_1` is the first order Bessel function of the first
    kind, :math: `r` is the radial distance from the maximum of the
    function (:math: `r = \\sqrt{(x - x_0)^2 + (y - y_0)^2}`), :math:`R`
    is the input ``radius`` parameter, and :math:`R_z =
    1.2196698912665045`.

    For an optical system, the radius of the first zero represents the
    limiting angular resolution and is approximately 1.22 * lambda / D,
    where lambda is the wavelength of the light and D is the diameter of
    the aperture.
    """

    amplitude = Parameter(default=1)
    x_0 = Parameter(default=0)
    y_0 = Parameter(default=0)
    radius = Parameter(default=1)
    eps = Parameter(default=0)
    _j1 = None

    def __init__(self, amplitude=amplitude.default, x_0=x_0.default,
                 y_0=y_0.default, radius=radius.default, eps=eps.default,
                 **kwargs):
        if self._j1 is None:
            try:
                from scipy.special import j1, jn_zeros
                self.__class__._j1 = j1
                self.__class__._rz = jn_zeros(1, 1)[0] / np.pi
            #add a ValueError here for python3 + scipy < 0.12
            except ImportError:
                raise ImportError("AiryDiskDiff2D model requires scipy > 0.11.")

        super(AiryDiskDiff2D, self).__init__(amplitude=amplitude,
                                             x_0=x_0,
                                             y_0=y_0,
                                             radius=radius,
                                             eps=eps,
                                             **kwargs)

    # Comment and methods copied from astropy v1.1.1
    # TODO: Why does this particular model have its own special __deepcopy__
    # and __copy__?  If it has anything to do with the use of the j_1 function
    # that should be reworked.
    def __deepcopy__(self, memo):
        new_model = self.__class__(self.amplitude.value, self.x_0.value,
                                   self.y_0.value, self.radius.value)
        return new_model

    def __copy__(self):
        new_model = self.__class__(self.amplitude.value, self.x_0.value,
                                   self.y_0.value, self.radius.value)
        return new_model

    @classmethod
    def evaluate(cls, x, y, amplitude, x_0, y_0, radius, eps):
        """Two dimensional Airy difference model function"""

        r = np.sqrt((x - x_0)**2 + (y - y_0)**2) / (radius / cls._rz)
        # Since r can be zero, we have to take care to treat that case
        # separately so as not to raise a numpy warning
        z = np.ones(r.shape)
        rt = np.pi * r[r > 0]
        z[r > 0] = (2.0 * cls._j1(rt) / rt -
                    2.0 * eps * cls._j1(eps * rt) / rt)**2
        z *= amplitude
        return z
