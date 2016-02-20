###############################################################################
# PlaneEffect
#
# DESCRIPTION
# The PlaneEffect functions are used to simulate effects that occur on all spectral
# layers equally, for example, sky rotation, telescope jitter, distortion, etc.
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
from copy import deepcopy

import numpy as np
import scipy.ndimage as spi

def gaussian_dist(x, mu, sig):
    p = np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    return p / np.sum(p)

def linear_dist(x):
    p = np.array([1.]*len(x))
    return p / np.sum(p)

def line_blur(arr, shift, kernel="gaussian", angle=0, pix_res=0.004):
    """
    Introduce a linear blur due to tracking error.

    Parameters
    ==========
    - arr: [2D array] the image
    - shift: [arcsec] how far in angular distance that the image has moved

    Optional parameters
    ===================
    - kernel: 'gaussian' - shift is the FWHM of the blur, approximating a random
                           walk in tracking error
              'linear' - shift is the length of the tracking blur with all
                         positions weighted equally, approximating no tracking
    - angle: [deg] the angle between image up and the zenith
    - pix_res: [arcsec] angular resolution of the pixels
    """

    # sample the shift at least every half pixel
    n = max(int(2 * shift) + 1, 3)
    if kernel == "gaussian":
        dr = np.linspace(-3 * shift, 3 * shift, 6 * n)
        weight = gaussian_dist(dr, 0, shift)
    else:
        dr = np.linspace(0, shift, n)
        weight = linear_dist(dr)

    dx = np.cos(np.deg2rad(ang)) * dr
    dy = np.sin(np.deg2rad(ang)) * dr

    tmp_arr = np.zeros(arr.shape)
    for x,y,w in zip(dx,dy,weight):
        tmp_arr += spi.shift(im, (dx, dy), order=1) * w

    return tmp_arr

def rotate_blur(arr, angle, kernel="gaussian"):
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

    angle_at_outer_pixel = np.rad2deg(np.arctan2(1, arr.shape[0] // 2))
    n = max(3, int(angle / angle_at_outer_pixel) + 1)

    if kernel == "gaussian":
        d_ang = np.linspace(-3 * angle, 3 * angle, max(2, 6*n))
        weight = gaussian(d_ang, 0, angle)
    else:
        d_ang = np.linspace(0, angle, n)
        weight = linear(d_ang)

    tmp_arr = np.zeros((q,q))
    for ang, w in zip(d_ang,weight):
        tmp_arr += spi.rotate(arr, ang, order=1, reshape=False) * w

    return tmp_arr


# def derotator(arr, angle, pix_res, kernel="gaussian"):


class CoordEffect(object):
    """
    """

    def __init__(self, dx, dy):
        self.dx = dx
        self.dy = dy

    def apply(self, x, y):
        return x + self.dx, y + self.dy



class ADC_Effect(CoordEffect):

    def __init__(lam, angle = 0, **kwargs):

        shift = atmospheric_refraction(lam, **kwargs)
        dx = shift * np.cos(np.deg2rad(angle))
        dy = shift * np.cos(np.deg2rad(angle))

        super(ADC_Effect, self).__init__(dx, dy)

class ArrayEffect(object):

    def __init__(self, x, y, weight):
        self.x = np.zeros()


    def apply(self, array):
        pass
