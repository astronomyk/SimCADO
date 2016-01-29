###############################################################################
# PSFCube
#
# Deliverables
#
#
#
#
#
#
#
# Classes:
#   
#
# Methods:
#   
#
#
#
#
#



# There are two types of psf object here:
#     - a cube
#     - a single psf image
# The cube is essentially a list of psf images, homogenized in size
# Should we have separate classes for these?
#
# Both PSF and PSFCube can be created from a single model or through
# convolution of a list of PSF components

from astropy import units as u
from astropy.convolution import (Gaussian2DKernel, AiryDisk2DKernel,
                                 Moffat2DKernel)
from astropy.convolution import convolve, convolve_fft
import utils


## These classes and functions are exported to the package
__all__ = ["PSFCube", "PSF"]

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
        z0 = float(config_dict['ATMO_ZENITH_DIST'])
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

        ## Rotate by the parallactice angle
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
            
        
class PSF(object):
    """Point spread function (single layer)

    Initialization methods:
    ----------------------

    PSF.gen_from_list : generates a PSF through convolution of
                        a list of psfs
    PSF.gen_analytic : generate PSF image from an analytic model 
    """

    def __init__(self):
        self.info = dict([])
        self.info['created'] = 'yes'

    @classmethod
    def gen_from_list(self, psf_list):
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
        
        
    @classmethod
    def gen_analytic(self, lam, kernel='airy', fwhm=0., res=1.*u.mas,
                     position=(0,0), padding=5, min_size=11, size=None):
        ## CHECK: How does this relate to the cube?
        """Generate a psf from an analytic model

        Parameters:
        ==========
        - lam : array of wavelength values
        - kernel : currently implemented are "airy" and "gauss"
                   For any other value, a delta peak is created.
        - fwhm : full-width at half maximum of the model. If 0, the theoretical
                 diffraction fwhm of the E-ELT is computed
        - res : ???
        - position : tuple giving the position of the delta peak 
                     (not for airy and gauss)
        - padding : increase output array by this margin outside of delta
                    peak
        - min_size : minimum size of output array (for airy and gauss)
        """
        
        self = PSF()
        self.lam = utils.unify(lam, u.um)
        self.kernel = kernel
        self.fwhm = utils.unify(fwhm, u.mas)
        self.res =  utils.unify(res, u.mas)

        ## convert sigma (gauss) to first zero (airy)
        gauss2airy = 2.76064  

        ## if lam does not have units, assume microns
        if self.fwhm.value == 0:
            ## TODO: telescope parameters as variables!
            self.fwhm = (1.22 * lam / (39.3 * u.m)).cgs * u.rad * (
                206264806 * u.mas / u.rad)

        ## create airy, gaussian or point source PSFs
        if "airy" in kernel:
            print "Making airy"
            self.info['description'] = "Airy PSF, FWHM = %.1f %s" \
                                       % (self.fwhm.value, self.fwhm.unit)
            length = int(max(16 * self.fwhm / self.res, min_size))
            if length % 2 == 0:
                length += 1
            n = (self.fwhm / 2.35) / self.res * gauss2airy
            self.psf = AiryDisk2DKernel(n, x_size=length, y_size=length).array
        elif "gauss" in kernel:
            print "Making gauss"
            self.info['description'] = "Gaussian PSF, FWHM = %.1f %s" \
                                       % (self.fwhm.value, self.fwhm.unit)
            length = int(max(8 * self.fwhm / self.res, min_size))
            if length % 2 == 0:
                length += 1

            n = (self.fwhm / 2.35) / self.res
            self.psf = Gaussian2DKernel(n, x_size=length, y_size=length,
                                        mode='oversample').array
        elif "moffat" in kernel:
            print "Making Moffat"
            ## CHECK: make this a parameter?
            beta = 4.765 ### Trujillo et al. 2001
            alpha = self.fwhm/(2 * np.sqrt(2**(1/beta) - 1))
            self.info['description'] = "Moffat PSF, FWHM = %.1f, alpha = %.1f"\
                                       % (self.fwhm.value, alpha.value)
            n = alpha/self.res
            if size is None:
                length = int(max(8 * self.fwhm / self.res, min_size))
            else:
                length = size
            if length % 2 == 0:
                length += 1
            if length > 100:
                mode = "linear_interp"
            else:
                mode = "oversample"
            self.psf = Moffat2DKernel(alpha.value, beta, x_size=length,
                                      y_size=length, mode=mode).array
        else:
            ## CHECK: What does this do?
            self.info['description'] = "Delta PSF, centred at (%.1f, %.1f)" \
                                       % position

            self.size = int(np.max(np.abs(position))) * 2 + padding
            if self.size % 2 == 0:
                self.size += 1
            self.x = self.size / 2 + position[0]
            self.y = self.size / 2 + position[1]

            x2 = self.x - int(self.x)
            x1 = 1. - x2
            y2 = self.y - int(self.y)
            y1 = 1. - y2

            n = np.zeros((self.size, self.size))
            n[int(self.y):int(self.y)+2, int(self.x):int(self.x)+2] = \
                np.array([[x1 * y1, x2 * y1], [x1 * y2, x2 * y2]])
            self.psf = n

        self.psf[self.psf <=0] = 1e-15
        self.psf = self.psf / np.sum(self.psf)
        return self

    
    def __repr__(self):
        return self.info['description']

    def resize(self, new_size):
        """Make all slices the same size

        The target shape is new_size x new_size
        """
        # make sure the new size is always an odd number
        if new_size % 2 == 0:
            new_size += 1

        for slice in self.cube:
            # create an empty square array of new size
            # get side length of the original array
            new_arr = np.zeros((new_size, new_size))
            aw, bw = slice.psf.shape[0], new_size

            # place old array in the centre of the empty array
            new_arr[bw/2 - aw/2 : bw/2 + aw/2 + 1, bw/2-aw/2 : bw/2 + aw/2 + 1] \
                = slice.psf
            slice.psf = new_arr

        self.size = [slice.psf.shape[0] for slice in self.cube]
            

## This is essentially a copy from Kieran recoded as a subclass
## Can we dispense with this class?
## Main benefit of having an extra class is automatic storage of
## atmospheric parameters.
#class PSF_ADC(PSF_Cube):
#    """Atmospheric dispersion stored as a "PSF" cube
#    """
#    
#    def __init__(self, lam_arr, res=1*u.mas, effectiveness=100,
#                 z0=60, temp=0, rel_hum=60, pres=750, nadir_angle=0,
#                 lat=-24.5, h=3064):
#        pass
                 
