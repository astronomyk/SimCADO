"""Atmosphere model"""
###############################################################################
# AtmosphereModel
#
# DESCRIPTION
# AtmosphereModel holds the transmission and emission curves as well as info
# on the environment variables: Temp, Pres, Rel Hum, z0. It should be able to
# decide whether a skycalc file is being read, of a series of 2 ascii files.
#
# This script also contains a series of helpful functions.
# - atmospheric_refraction() computes the angular difference between the
#     apparent position of a star seen from the ground and its true position.
#
# Classes:
#   AtmosphereModel(object)
#   - __init__(UserCommands)
#   - reshape(lam)
#   -
#
# Methods:
#   atmospheric_refraction(lam, z0=60, temp=0, rel_hum=60, pres=750, lat=-24.5,
#                          h=3064):
#
#
#
#
#


import numpy as np
from astropy import units as u
from astropy.io import fits
import SimCADO.SpectralCurve as sc

## These classes and functions are exported to the package
__all__ = ["AtmosphereModel", "atmospheric_refraction"]

class AtmosphereModel(object):
    """Class holding information about the atmosphere"""

    def __init__(self, fname, res=0.001*u.um):
        ## CHECK: What's the purpose of res?
        ## TODO: Identify fname from input parameters?
        self.info = dict([])

        self.info["fname"] = fname
        self.info["res"] = res

        hdulist = fits.open(fname)
        data = hdulist[1].data
        hdulist.close()

        ## CHECK: Do these have to go into the top level?
        lam = data["lam"]
        flux = data["flux"]
        trans = data["trans"]
        val_unit = u.Unit("ph/s/m2/micron/arcsec2")

        self.emission = sc.EmissionCurve(lam=lam, val=flux,
                                         lam_unit=u.um, val_unit=val_unit,
                                         Type="Atmospheric emission")
        self.transmission = sc.TransmissionCurve(lam=lam, val=trans,
                                                 lam_unit=u.um, val_unit=1.,
                                                 Type="Atmospheric transmission")



    def __repr__(self):
        return "Ich bin ein Atmosphere_Model:\n" + str(self.info)


def atmospheric_refraction(lam, z0=60, temp=0, rel_hum=60, pres=750,
                           lat=-24.5, h=3064):
    """Compute atmospheric refraction

    The function computes the angular difference between the apparent position
    of a star seen from the ground and its true position.

    Parameters:
    ==========
    - lam : wavelength [micron]
    - z0 : zenith distance [degrees]
    - temp : ground temperature [deg C]
    - rel_hum : relative humidity [%]
    - pres : air pressure [millibar]
    - lat : latitude [degrees] - set for Cerro Armazones: 24.5 deg South
    - h : height above sea level [meters] - 3064 m

    Output:
    ======
    - ang : angle between real position and refracted position [arcsec]

    References:
    ==========
    See Stone 1996 and the review by S. Pedraz -
    http://www.caha.es/newsletter/news03b/pedraz/newslet.html
    """

    #### need T, P, RH for Ps, Pw Pa
    T = 273.15 + temp

    Ps = -10474. + 116.43 * T - 0.43284 * T**2 + 0.0005384 * T**3
    Pw = Ps * rel_hum / 100.
    Pa = pres - Pw

    #### need n0 for gamma
    sig = 1. / lam
    Da = (Pa / T) * (1. + Pa * (57.9E-8 - 0.0009325 / T + 0.25844 / T**2))
    Dw = (Pw / T) * (1. + Pw * (1. + 3.7E-4 * Pw) *
                     (-2.37321E-3 + 2.23366 / T - 710.792 / T**2
                     + 77514.1 / T**3))

    na = Da * (2371.34 + 683939.7 / (130. - sig**2) + 4547.3 / (38.9 - sig**2))
    nw = Dw * (6487.31 + 58.058 * sig**2 - 0.7115 * sig**4 + 0.08851 * sig**6)
    n0 = 1E-8 * (na + nw) + 1.

    #### need gamma, kappa and beta for R
    g = n0 - 1.
    b = 0.001254 * (273.15 + temp) / 273.15
    #k = 1
    k = 1. + 0.005302 * np.sin(lat / 57.29578)**2 \
        - 0.00000583 * np.sin(2. * lat / 57.29578)**2 - 0.000000315 * h

    R = k * g * (1 - b) * np.tan(z0 / 57.29578) \
        - k * g * (b - g/2.) * np.tan(z0/57.29578)**3

    # the refraction is the side of a triangle, although the triangle
    # side makes an arc across the sky.
    # We want the angle that this triangle side is subtending
    # Using the small angle approximation (which is in radians),
    # we can get the angle of refraction

    ang = R * 3600 * 57.29578

    return ang
