'''Classes to describe PSFs in SimCADO'''
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
from astropy.convolution import convolve_fft

import numpy as np

import SimCADO.utils as utils
from SimCADO.atmosphere_model import atmospheric_refraction

## These classes and functions are exported to the package
__all__ = ["PSFCube", "PSF"]

class PSFCube(object):
    """Class holding wavelength dependent point spread function"""

    def __init__(self):
        self.info = dict([])
        self.info['created'] = 'yes'

        # Attributes that are filled in by methods
        # TODO: Put some into self.info
        self.kernel = None
        self.cube = None
        self.lam = None
        self.res = None
        self.psf = None
        self.padding = None
        self.lam_arr = None
        self.size_orig = None
        self.fwhm = None

    def __repr__(self):
        return self.info['description']

    @classmethod
    def gen_adc(cls, lam_arr, config_dict, res=1*u.mas, nadir_angle=0):
        # TODO: Reformulate doc string
        '''Generate a PSF cube to describe atmospheric dispersion

Each layer in the cube has a delta peak psf at a position corresponding
to the relative atmospheric refraction of that wavelength.
        '''
        # TODO: Remove nadir angle (=parallactic angle?)
        # TODO: add info
        effectiveness = float(config_dict['INST_ADC_EFFICIENCY'])/100. ###
        lat = float(config_dict['SCOPE_LATITUDE'])
        alt = float(config_dict['SCOPE_ALTITUDE'])
        rel_hum = float(config_dict['ATMO_REL_HUMIDITY'])

        ## CHECK: better called OBS_ZENITH_DIST
        zendist = float(config_dict['ATMO_ZENITH_DIST'])
        temp = float(config_dict['ATMO_TEMPERATURE'])
        pres = float(config_dict['ATMO_PRESSURE'])

        lam_arr = utils.unify(lam_arr, u.um)
        nadir_angle = utils.unify(nadir_angle, u.deg)
        res = utils.unify(res, u.mas)

        ## get the angle shift for each slice
        angle_shift = [atmospheric_refraction(lam, zendist, temp, rel_hum,
                                              pres, lat, alt)
                       for lam in lam_arr.value]
        angle_shift = angle_shift * u.arcsec

        ## convert angle shift into number of pixels
        ## pixel shifts are defined with respect to last slice
        pixel_shift = (angle_shift - angle_shift[-1]).to(u.mas) / res

        ## Rotate by the parallactice angle
        xpos = -pixel_shift * np.sin(nadir_angle.to(u.rad)) \
                 * (1. - effectiveness)
        ypos = -pixel_shift * np.cos(nadir_angle.to(u.rad)) \
                 * (1. - effectiveness)

        adc_cube = cls.gen_cube(lam_arr, kernel='adc',
                                position=(xpos, ypos),
                                res=res)
        return adc_cube

    @classmethod
    def gen_cube(cls, lam_arr, fwhm=0, res=1.*u.mas, kernel='gauss',
                 position=(0, 0), padding=5, min_size=71):
        """Generate a cube full of psf_gen objects

        An astropy unit array is required for the central wavelength
        of each slice.
        The FWHM can be specified as a blanket value, or an array of
        FWHMs for each slice
        """

        ## Get everything in mas for the spatial values and um for the
        ## spectral values
        psfcube = cls()

        if 'airy' in kernel:
            psfcube.info['description'] = "PSF cube, Airy"
        elif 'gauss' in kernel:
            psfcube.info['description'] = "PSF cube, Gauss"
        else:
            psfcube.info['description'] = "PSF cube, delta peak"

        psfcube.lam_arr = utils.unify(lam_arr, u.um)

        # if fwhm is scalar, generate a constant array
        psfcube.fwhm = utils.unify(fwhm, u.mas, len(lam_arr))
        psfcube.res = utils.unify(res, u.mas)
        psfcube.kernel = kernel
        psfcube.padding = padding

        ## Create a cube along the spectral dimension. Call psf_gen for each
        ## lambda value and corresponding fwhm
        psfcube.cube = []

        ## if only one position, generate constant arrays for x and y
        xpos, ypos = position[0], position[1]
        if not hasattr(xpos, "__len__"):
            xpos = xpos * np.ones(len(lam_arr))
        if not hasattr(ypos, "__len__"):
            ypos = ypos * np.ones(len(lam_arr))

        for i in range(len(psfcube.lam_arr)):
            psfcube.cube += [PSF.gen_analytic(psfcube.lam_arr[i],
                                              fwhm=psfcube.fwhm[i],
                                              res=psfcube.res,
                                              kernel=psfcube.kernel,
                                              position=(xpos[i], ypos[i]),
                                              padding=padding,
                                              min_size=min_size)]

        psfcube.fwhm = [i.fwhm for i in psfcube.cube]
        psfcube.size_orig = [i.psf.shape[0] for i in psfcube.cube]

        ## CHECK: Do we have to resize the slices to a common shape?
        ## Comment from Kieran's original code
        # resample the slices so that they are all the same size
        # psfcube.resample(np.max(psfcube.size_orig))

	# !!!!!!!!!!! It seems that we don't need the resampling,
        # so long as all arrays are odd numbers in length !!!!!!!!!!!!!!

        ## TODO: Add some more cube info

        return psfcube


    @classmethod
    def gen_from_list(cls, psf_list):
        """Generate a master psf cube through convolution of a list of psfs"""

        if not hasattr(psf_list, "__len__"):
            psf_list = [psf_list]

        psfcube = cls()

        psfcube.info['description'] = "Master psf cube from list"
        for i in range(len(psf_list)):
            psfcube.info['PSF%02d' % (i+1)] = psf_list[i].info['description']

        ## Check that the wavelengths are equal
        lam_list = [psf.lam_arr for psf in psf_list]
        if not all([all(lam == lam_list[0]) for lam in lam_list[1:]]):
            raise ValueError("Wavelength arrays of psf cubes are not equal")
        psfcube.lam = lam_list[0]

        ## Check that the resolutions are equal
        res_list = [psf.res for psf in psf_list]
        if not all([res == res_list[0] for res in res_list[1:]]):
            raise ValueError("Resolutions of psf cubes are not equal")
        psfcube.res = res_list[0]

        ## Combine each layer
        # CHECK: What's this padding for? size is not defined
        ##padding = np.sum(np.max([psf.size for psf in psf_list], axis=1))

        psfcube.psf = PSFCube.gen_cube(psfcube.lam, res=psfcube.res,
                                       kernel='dot')
        #, min_size=min_size, padding=padding)

        psfcube.cube = range(len(psfcube.lam))
        for i in range(len(psfcube.lam)):
            psfcube.cube[i] = PSF.gen_from_list([psf.cube[i]
                                                 for psf in psf_list])

        return psfcube


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

        # Attributes that are assigned by methods
        self.lam = None
        self.kernel = None
        self.res = None
        self.psf = None
        self.x = None
        self.y = None
        self.fwhm = None
        self.size = None

    @classmethod
    def gen_from_list(cls, psf_list):
        """Generate a master psf through convolution of a list of psfs"""
        import copy

        if not hasattr(psf_list, "__len__"):
            psf_list = [psf_list]

        psfcube = copy.deepcopy(psf_list[0])
        psfcube.info['description'] = "Master psf from list"
        psfcube.info['PSF01'] = psf_list[0].info['description']
        for i, psf in enumerate(psf_list[1:]):
            psfcube.psf = convolve_fft(psfcube.psf, psf.psf)
            psfcube.info['PSF%02d' % (i+2)] = psf.info['description']

        return psfcube


    @classmethod
    def gen_analytic(cls, lam, kernel='airy', fwhm=0., res=1.*u.mas,
                     position=(0, 0), padding=5, min_size=11, size=None):
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

        psf = PSF()
        psf.lam = utils.unify(lam, u.um)
        psf.kernel = kernel
        psf.fwhm = utils.unify(fwhm, u.mas)
        psf.res = utils.unify(res, u.mas)

        ## convert sigma (gauss) to first zero (airy)
        gauss2airy = 2.76064

        ## if lam does not have units, assume microns
        if psf.fwhm.value == 0:
            ## TODO: telescope parameters as variables!
            psf.fwhm = (1.22 * lam / (39.3 * u.m)).cgs * u.rad * (
                206264806 * u.mas / u.rad)

        ## create airy, gaussian or point source PSFs
        if "airy" in kernel:
            print "Making airy"
            psf.info['description'] = "Airy PSF, FWHM = %.1f %s" \
                                       % (psf.fwhm.value, psf.fwhm.unit)
            length = int(max(16 * psf.fwhm / psf.res, min_size))
            if length % 2 == 0:
                length += 1
            n = (psf.fwhm / 2.35) / psf.res * gauss2airy
            psf.psf = AiryDisk2DKernel(n, x_size=length, y_size=length).array
        elif "gauss" in kernel:
            print "Making gauss"
            psf.info['description'] = "Gaussian PSF, FWHM = %.1f %s" \
                                       % (psf.fwhm.value, psf.fwhm.unit)
            length = int(max(8 * psf.fwhm / psf.res, min_size))
            if length % 2 == 0:
                length += 1

            n = (psf.fwhm / 2.35) / psf.res
            psf.psf = Gaussian2DKernel(n, x_size=length, y_size=length,
                                       mode='oversample').array
        elif "moffat" in kernel:
            print "Making Moffat"
            ## CHECK: make this a parameter?
            beta = 4.765 ### Trujillo et al. 2001
            alpha = psf.fwhm/(2 * np.sqrt(2**(1/beta) - 1))
            psf.info['description'] = "Moffat PSF, FWHM = %.1f, alpha = %.1f"\
                                       % (psf.fwhm.value, alpha.value)
            n = alpha/psf.res
            if size is None:
                length = int(max(8 * psf.fwhm / psf.res, min_size))
            else:
                length = size
            if length % 2 == 0:
                length += 1
            if length > 100:
                mode = "linear_interp"
            else:
                mode = "oversample"
            psf.psf = Moffat2DKernel(alpha.value, beta, x_size=length,
                                     y_size=length, mode=mode).array
        else:
            ## CHECK: What does this do?
            psf.info['description'] = "Delta PSF, centred at (%.1f, %.1f)" \
                                       % position

            psf.size = int(np.max(np.abs(position))) * 2 + padding
            if psf.size % 2 == 0:
                psf.size += 1
            psf.x = psf.size / 2 + position[0]
            psf.y = psf.size / 2 + position[1]

            x2 = psf.x - int(psf.x)
            x1 = 1. - x2
            y2 = psf.y - int(psf.y)
            y1 = 1. - y2

            n = np.zeros((psf.size, psf.size))
            n[int(psf.y):int(psf.y)+2, int(psf.x):int(psf.x)+2] = \
                np.array([[x1 * y1, x2 * y1], [x1 * y2, x2 * y2]])
            psf.psf = n

        psf.psf[psf.psf <= 0] = 1e-15
        psf.psf = psf.psf / np.sum(psf.psf)
        return psf


    def __repr__(self):
        return self.info['description']

    def resize(self, new_size):
        """Make all slices the same size

        The target shape is new_size x new_size
        """
        # make sure the new size is always an odd number
        if new_size % 2 == 0:
            new_size += 1

        for psf in self.cube:
            # create an empty square array of new size
            # get side length of the original array
            new_arr = np.zeros((new_size, new_size))
            aw, bw = psf.psf.shape[0], new_size

            # place old array in the centre of the empty array
            new_arr[bw/2 - aw/2 : bw/2 + aw/2 + 1,
                    bw/2-aw/2 : bw/2 + aw/2 + 1] = psf.psf
            psf.psf = new_arr

        ## CHECK: pylint says PSF has no member 'cube'
        self.size = [psf.psf.shape[0] for psf in self.cube]


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
