"""
Helper functions for SimCADO
"""
###############################################################################
# utils.py
#
# DESCRIPTION
#
#
# Classes:
#
#
# Methods:
#  read_config(config_file)
#  update_config(config_file, config_dict)
#  unify(x, unit, length=1)
#  parallactic_angle(ha, de, lat=-24.589167)
#  parallactic_angle_2(ha, de, lat=-24.589167)
#  moffat(r, alpha, beta)
#


import numpy as np
from astropy import units as u
from astropy.io import fits

## These functions are exported to the package
__all__ = ["read_config", "update_config", "unify", "parallactic_angle",
           "poissonify", "atmospheric_refraction", "nearest", "add_keyword"]


def msg(cmds, message, level=3):
    """
    Prints a message based on the level of verbosity given in cmds

    Parameters
    ----------
    cmds : UserCommands
        just for the SIM_VERBOSE and SIM_MESSAGE_LEVEL keywords
    message : str
        message to be printed
    level : int, optional
        all messages with level <= SIM_MESSAGE_LEVEL are printed. I.e. level=5
        messages are not important, level=1 are very important
    """
    if cmds["SIM_VERBOSE"] == "yes" and level <= cmds["SIM_MESSAGE_LEVEL"]:
        print(message)


## CHECK: Turn config into a class? (subclass of dict, if anything)
def read_config(config_file):
    """
    Read in a SimCADO configuration file

    A configuration file in the SExtractor format:
       'PARAMETER    Value    # Comment'

    Parameters
    ----------
    config_file : str
        the filename of the .config file

    Returns
    -------
    config_dict : dict
        A dictionary with keys 'PARAMETER' and values 'Value'.

    Notes
    -----
    The values of the dictionary are strings and will have to be converted to
    the appropriate data type as they are needed.
    """

    import re

    # Read the file into a list of strings
    config_dict = dict()

    # remove lines that are all spaces or spaces + '#'
    # these are the regular expressions
    isempty = re.compile(r'^\s*$')
    iscomment = re.compile(r'^\s*#')

    with open(config_file, 'r') as fp1:
        for line in fp1:
            if isempty.match(line):
                continue
            if iscomment.match(line):
                continue

            line = line.rstrip()             # remove trailing \n
            content = line.split('#', 1)[0]  # remove comment
            param, value = content.split(None, 1)

            # Convert to number if possible
            try:
                config_dict[param] = float(value.strip())
            except ValueError:
                config_dict[param] = value.strip()

            # Convert string "none" to python None
            if isinstance(value, str) and value.lower() == "none":
                config_dict[param] = None

    return config_dict


def update_config(config_file, config_dict):
    """
    Update a SimCADO configuration dictionary

    A configuration file in the SExtractor format:
         'PARAMETER    Value    # Comment'
    an existing configuration dictionary.

    Parameters
    ----------
    config_file : str
        the filename of the .config file

    Returns
    -------
    config_dict : dict
        A dictionary with keys 'PARAMETER' and values 'Value'.

    Returns:
    -------
    config_dict : dict
        A dictionary with keys 'PARAMETER' and values 'Value'.

    Notes
    -----
    the values of the dictionary are strings and will have
    to be converted to the appropriate data type as they are needed.
    """
    config_dict.update(read_config(config_file))

    return config_dict

def unify(x, unit, length=1):
    """
    Convert all types of input to an astropy array/unit pair

    Parameters
    ----------
    x : int, float, np.ndarray, astropy.Quantity
        The array to be turned into an astropy.Quantity
    unit : astropy.Quantity
        The units to attach to the array
    length : int, optional
        If `x` is a scalar, and the desired output is an array with `length`

    Returns
    -------
    y : astropy.Quantity
    """

    print(type(x))

    if isinstance(x, u.quantity.Quantity):
        if isinstance(x.value, np.ndarray):
            y = x.to(unit)
        elif length == 1:
            y = x.to(unit)
        else:
            y = ([x.value] * length * x.unit).to(unit)
    else:
        if hasattr(x, "__len__"):
            y = x * unit
        elif length == 1:
            y = x * unit
        else:
            y = [x] * length * unit

    return y


### Just a sketch
def parallactic_angle(ha, de, lat=-24.589167):
    """
    Compute the parallactic angle

    Parameters
    ----------
    ha, de : float
        hour angle and declination of observation point
    lat : float
        latitude of observatory

    Returns
    -------
    parang : float
        The parallactic angle
    """

    # Convert observation point, pole and zenith to cartesian coordinates
    x_pole = [0, 0, -1]
    print(x_pole)

    poledist = (90. - lat)/180 * np.pi  # from north pole
    x_zenith = [0, np.sin(poledist), np.cos(poledist)]
    print(x_zenith)

    obsdist = (90. - de)/180 * np.pi
    ha = ha/180. * np.pi
    x_obs = [np.sin(obsdist) * np.sin(ha), np.sin(obsdist) * np.cos(ha),
             np.cos(obsdist)]
    print(x_obs)

    # normals to the great circles
    n_pole = np.cross(x_obs, x_pole)
    n_pole = n_pole/np.sqrt((n_pole**2).sum())
    n_zenith = np.cross(x_obs, x_zenith)
    n_zenith = n_zenith/np.sqrt((n_zenith**2).sum())

    # Angle between the great circles
    ## CHECK: exact definition, modulo 180 deg?
    parang = np.arccos((n_pole * n_zenith).sum()) * 180. / np.pi

    return parang

# From Filippenko1982
def parallactic_angle_2(ha, de, lat=-24.589167):
    """
    Compute the parallactic angle

    Parameters
    ----------
    ha, de : float
        hour angle and declination of observation point
    lat : float
        latitude of observatory

    Returns
    -------
    parang : float
        The parallactic angle

    References
    ----------
    Filippenko (1982)
    """

    hadeg = ha
    ha = ha/180. * np.pi
    de = de/180. * np.pi
    lat = lat/180. * np.pi

    sineta = np.sin(ha) * np.cos(lat) / \
             np.sqrt(1. - (np.sin(lat) * np.sin(de) +
                           np.cos(lat) * np.cos(de) * np.cos(ha))**2)

    eta = np.arcsin(sineta) * 180./np.pi
    if hadeg >= 0 and hadeg <= 45:
        eta = 180. - eta
    elif hadeg < 0 and hadeg >= -45:
        eta = 180. + eta
    elif hadeg < -45:
        eta = - eta
    return eta


def moffat(r, alpha, beta):
    """
    !!Unfinished!! Return a Moffat function

    Parameters
    ----------
    r
    alpha
    beta

    Returns
    -------
    eta
    """
    return (beta - 1)/(np.pi * alpha**2) * (1 + (r/alpha)**2)**(-beta)


def poissonify(arr):
    """
    Add a realisation of the poisson process to the array 'arr'.

    Parameters
    ----------
    arr : np.ndarray
        The input array which needs a Poisson distribution applied to items

    Returns
    -------
    arr : np.ndarray
        The input array, but with every pixel altered according to a poisson
        distribution
    """
    return np.random.poisson(arr).astype(np.float32)


def atmospheric_refraction(lam, z0=60, temp=0, rel_hum=60, pres=750,
                           lat=-24.5, h=3064):
    """Compute atmospheric refraction

    The function computes the angular difference between the apparent position
    of a star seen from the ground and its true position.

    Parameters
    ----------
    lam : float, np.ndarray
        [um] wavelength bin centres
    z0 : float, optional
        [deg] zenith distance. Default is 60 deg from zenith
    temp : float, optional
        [deg C] ground temperature. Default is 0 deg C
    rel_hum : float, optional
        [%] relative humidity. Default is 60%
    pres : float, optional
        [millibar] air pressure. Default is 750 mbar
    lat : float, optional
        [deg] latitude. Default set for Cerro Armazones: 24.5 deg South
    h : float, optional
        [m] height above sea level. Default is 3064 m

    Returns
    -------
    ang : float, np.ndarray
        [arcsec] angle between real position and refracted position

    References
    ----------
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
    k = 1. + 0.005302 * np.sin(np.deg2rad(lat))**2 \
        - 0.00000583 * np.sin(2. * np.deg2rad(lat))**2 - 0.000000315 * h

    R = k * g * (1 - b) * np.tan(np.deg2rad(z0)) \
        - k * g * (b - g/2.) * np.tan(np.deg2rad(z0))**3

    # the refraction is the side of a triangle, although the triangle
    # side makes an arc across the sky.
    # We want the angle that this triangle side is subtending
    # Using the small angle approximation (which is in radians),
    # we can get the angle of refraction

    ang = np.rad2deg(R * 3600)

    return ang


def nearest(arr, val):
    """
    Return the index of the value from 'arr' which is closest to 'val'

    Parameters
    ----------
    arr : np.ndarray, list, tuple
        Array to be searched
    val : float, int
        Value to find in `arr`

    Returns
    -------
    i : int
        index of array where the nearest value to `val` is
    """
    if isinstance(val, (list, tuple, np.ndarray)):
        arr = np.array(arr)
        return [nearest(arr, i) for i in val]

    d = arr - val
    i = np.where(abs(d) == np.min(abs(d)))[0][0]
    return i


def add_keyword(filename, keyword, value, comment="", ext=0):
    """
    Add a keyword, value pair to an extension header in a FITS file

    Parameters
    ----------
    filename : str
        Name of the FITS file to add the keyword to
    keyword : str
    value : str, float, int
    comment : str
    ext : int, optional
        The fits extension index where the keyword should be added. Default is 0
    """
    f = fits.open(filename, mode="update")
    f[ext].header[keyword] = (value, comment)
    f.flush()
    f.close()
