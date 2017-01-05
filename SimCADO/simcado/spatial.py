'''TODO: Insert docstring'''

###############################################################################
# PlaneEffect
#
# DESCRIPTION
# The PlaneEffect functions are used to simulate effects that occur on all
# spectral layers equally, for example, sky rotation, telescope jitter,
# distortion, etc.
# To do this in the most general way possible, a PlaneEffect function contains
# 3 planes representing the deviation in position from an ideal optical train.
# The values in each of the 3 planes represent the distance the pixel should
# move in the x and y directions, and a weighting value.
#
# Several functions generate the various effects that occur. For example:
# - Rotation
# - Distortion
# - Translation
# - FlatField
# Some PlaneEffects only need to act on the positions of the incoming photons,
# e.g. ADC, while others are applicable to the whole array, e.g. Distortion,
# Flat field. As each PlaneEffect
# - CoordEffect
# - ArrayEffect
#
# As each PlaneEffect
#
# Classes:
#  CoordEffect
#  ArrayEffect
#
# Subclasses:
#  Rotation(ArrayEffect)
#  Distortion(ArrayEffect)
#  FlatField(ArrayEffect)
#  Translation(CoordEffect)

#
# Methods:
#
#
#

#from copy import deepcopy   ## Not used (OC)

import os
from .utils import __pkg_dir__

import warnings, logging

import numpy as np
import scipy.ndimage as spi
from scipy.signal import fftconvolve

from astropy.convolution import convolve_fft, Gaussian2DKernel

try:
    import simcado.utils as utils
except ImportError:
    import utils


__all__ = ["tracking", "derotator", "wind_jitter", "adc_shift"]


def _gaussian_dist(x, mu, sig):
    p = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    return p / np.sum(p)

def _linear_dist(x):
    p = np.array([1.]*len(x))
    return p / np.sum(p)

def _line_blur(arr, shift, kernel="gaussian", angle=0):
    """
    Introduce a linear blur due to tracking error.

    Parameters
    ----------
    - arr: [2D array] the image
    - shift: [pixel] how many pixels the image has moved

    Optional parameters
    -------------------
    - kernel: 'gaussian' - shift is the FWHM of the blur, approximating a random
                           walk in tracking error
              'linear' - shift is the length of the tracking blur with all
                         positions weighted equally, approximating no tracking
    - angle: [deg] the angle between image up and the zenith
    """

    # sample the shift at least every half pixel
    n = max(int(2 * shift) + 1, 3)
    if kernel == "gaussian":
        dr = np.linspace(-3 * shift, 3 * shift, 6 * n)
        weight = _gaussian_dist(dr, 0, shift)
    else:
        dr = np.linspace(0, shift, n)
        weight = _linear_dist(dr)

    dx = np.cos(np.deg2rad(angle)) * dr
    dy = np.sin(np.deg2rad(angle)) * dr

    tmp_arr = np.zeros(arr.shape)
    for x, y, w in zip(dx, dy, weight):  # TODO: x, y unused? (OC)
        tmp_arr += spi.shift(arr, (x, y), order=1) * w

    return tmp_arr

def _rotate_blur(arr, angle, kernel="gaussian"):
    """
    Introduce a rotational blur due to derotator error.

    Parameters
    ==========
    - arr: [2D array] the image
    - angle: [deg] the angle or rotation

    Optional parameters
    ===================
    - kernel: 'gaussian' - angle is the FWHM of the blur, approximating a random
                           walk in tracking error
              'linear' - shift is the length of the tracking blur with all
                         positions weighted equally, approximating no tracking
    """

    ang_at_cnr_pix = np.rad2deg(np.arctan2(1, np.sqrt(2) * arr.shape[0] // 2))
    n = max(3, int(angle / ang_at_cnr_pix) + 1)

    if kernel == "gaussian":
        d_ang = np.linspace(-3 * angle, 3 * angle, max(2, 6*n))
        weight = _gaussian_dist(d_ang, 0, angle)
    else:
        d_ang = np.linspace(0, angle, n)
        weight = _linear_dist(d_ang)

    tmp_arr = np.zeros(arr.shape)
    for ang, w in zip(d_ang, weight):
        tmp_arr += spi.rotate(arr, ang, order=1, reshape=False) * w

    return tmp_arr


def tracking(arr, cmds):
    """
    A method to simulate tracking errors
    ===== Currently a place holder with minimum functionality =========
    !! TODO, work out the shift during the DIT for the object RA, DEC etc !!
    """
    if cmds["SCOPE_DRIFT_DISTANCE"] > 0.:
        pix_res = cmds["SIM_DETECTOR_PIX_SCALE"] / cmds["SIM_OVERSAMPLING"]
        kernel = cmds["SCOPE_DRIFT_PROFILE"]
        shift = cmds["SCOPE_DRIFT_DISTANCE"] / pix_res

        return _line_blur(arr, shift, kernel=kernel, angle=0)
    else:
        return arr


def derotator(arr, cmds):
    """
    A method to simulate field rotation in case the derotator is <100% effective
    ===== Currently a place holder with minimum functionality =========
    !! TODO, work out the rotation during the DIT for the object RA, DEC etc !!
    """
    if cmds["INST_DEROT_PERFORMANCE"] < 100.:
        eff = 1. - (cmds["INST_DEROT_PERFORMANCE"] / 100.)
        kernel = cmds["INST_DEROT_PROFILE"]
        angle = eff * cmds["OBS_EXPTIME"] * 15 / 3600.
        
        edge_smear = angle / cmds.pix_res
        if edge_smear > 50:
            print("The smear at the detector edge is large:", edge_smear, "[px]")
        
        return _rotate_blur(arr, angle, kernel=kernel)
    else:
        return arr


def wind_jitter(arr, cmds):
    """
    A method to simulate wind jitter
    ===== Currently a place holder with minimum functionality =========
    !! TODO, get the read spectrum for wind jitter !!
    !! Add in an angle parameter for the ellipse   !!
    """
    pix_res = cmds["SIM_DETECTOR_PIX_SCALE"] / cmds["SIM_OVERSAMPLING"]
    fwhm = cmds["SCOPE_JITTER_FWHM"] / pix_res
    n = (fwhm / 2.35)
    kernel = Gaussian2DKernel(n, mode="oversample")

    return fftconvolve(arr, kernel, mode="same")
    #return convolve_fft(arr, kernel, allow_huge=True)


def adc_shift(cmds):
    """Generates a list of x and y shifts from a commands object"""

    para_angle = cmds["OBS_PARALLACTIC_ANGLE"]
    effectiveness = cmds["INST_ADC_PERFORMANCE"] / 100.

    ## get the angle shift for each slice
    angle_shift = [utils.atmospheric_refraction(lam,
                                                cmds["OBS_ZENITH_DIST"],
                                                cmds["ATMO_TEMPERATURE"],
                                                cmds["ATMO_REL_HUMIDITY"],
                                                cmds["ATMO_PRESSURE"],
                                                cmds["SCOPE_LATITUDE"],
                                                cmds["SCOPE_ALTITUDE"])
                   for lam in cmds.lam_bin_centers]

    ## convert angle shift into number of pixels
    ## pixel shifts are defined with respect to last slice
    rel_shift = (angle_shift - angle_shift[-1])
    if np.max(np.abs(rel_shift)) > 1000:
        raise ValueError("Pixel shifts too great (>1000), check units")

    ## Rotate by the paralytic angle
    x = -rel_shift * np.sin(np.deg2rad(para_angle)) * (1. - effectiveness)
    y = -rel_shift * np.cos(np.deg2rad(para_angle)) * (1. - effectiveness)

    ## return values are in [arcsec]
    return x, y
