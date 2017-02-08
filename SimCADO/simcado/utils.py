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
# Functions:
#  read_config(config_file)
#  update_config(config_file, config_dict)
#  unify(x, unit, length=1)
#  parallactic_angle(ha, de, lat=-24.589167)
#  parallactic_angle_2(ha, de, lat=-24.589167)
#  moffat(r, alpha, beta)
#
#
import os
import inspect
from glob import glob

try:
    import wget
except ImportError:
    print("Package wget is not available. simcado.get_extras() will not work.")

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.io import ascii as ioascii


__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))

#__all__ = []
#__all__ = ["unify", "parallactic_angle", "poissonify",
#           "atmospheric_refraction", "nearest", "add_keyword"]


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
        If ``x`` is a scalar, and the desired output is an array with ``length``

    Returns
    -------
    y : astropy.Quantity
    """

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


# Just a sketch
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

    poledist = (90. - lat)/180 * np.pi  # from north pole
    x_zenith = [0, np.sin(poledist), np.cos(poledist)]

    obsdist = (90. - de)/180 * np.pi
    ha = ha/180. * np.pi
    x_obs = [np.sin(obsdist) * np.sin(ha), np.sin(obsdist) * np.cos(ha),
             np.cos(obsdist)]

    # normals to the great circles
    n_pole = np.cross(x_obs, x_pole)
    n_pole = n_pole/np.sqrt((n_pole**2).sum())
    n_zenith = np.cross(x_obs, x_zenith)
    n_zenith = n_zenith/np.sqrt((n_zenith**2).sum())

    # Angle between the great circles
    # CHECK: exact definition, modulo 180 deg?
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

    # need T, P, RH for Ps, Pw Pa
    T = 273.15 + temp

    Ps = -10474. + 116.43 * T - 0.43284 * T**2 + 0.0005384 * T**3
    Pw = Ps * rel_hum / 100.
    Pa = pres - Pw

    # need n0 for gamma
    sig = 1. / lam
    Da = (Pa / T) * (1. + Pa * (57.9E-8 - 0.0009325 / T + 0.25844 / T**2))
    Dw = (Pw / T) * (1. + Pw * (1. + 3.7E-4 * Pw) *
                     (-2.37321E-3 + 2.23366 / T - 710.792 / T**2 +
                      77514.1 / T**3))

    na = Da * (2371.34 + 683939.7 / (130. - sig**2) + 4547.3 / (38.9 - sig**2))
    nw = Dw * (6487.31 + 58.058 * sig**2 - 0.7115 * sig**4 + 0.08851 * sig**6)
    n0 = 1E-8 * (na + nw) + 1.

    # need gamma, kappa and beta for R
    g = n0 - 1.
    b = 0.001254 * (273.15 + temp) / 273.15
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

    # return value is in arcsec
    return ang


def nearest(arr, val):
    """
    Return the index of the value from 'arr' which is closest to 'val'

    Parameters
    ----------
    arr : np.ndarray, list, tuple
        Array to be searched
    val : float, int
        Value to find in ``arr``

    Returns
    -------
    i : int
        index of array where the nearest value to ``val`` is
    """
    if isinstance(val, (list, tuple, np.ndarray)):
        arr = np.array(arr)
        return [nearest(arr, i) for i in val]

    return np.argmin(abs(arr - val))


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
        The fits extension index where the keyword should be added.
        Default is 0
    """
    f = fits.open(filename, mode="update")
    f[ext].header[keyword] = (value, comment)
    f.flush()
    f.close()


############# Check the server for data extras
def download_file(url, save_dir=os.path.join(__pkg_dir__, "data")):
    """
    Download the extra data that aren't in the SimCADO package
    """

    local_filename = os.path.join(save_dir, url.split('/')[-1])
    try:
        temp_file = wget.download(url,
                                  out=wget.tempfile.mktemp(dir=save_dir, 
                                                           suffix='.tmp'),
                                  bar=wget.bar_adaptive)
        print("\n")
        if os.path.exists(local_filename): 
            os.remove(local_filename)
        os.rename(temp_file, local_filename)
    except wget.ulib.HTTPError:
        print(url + " not found")

    return local_filename


def get_extras():
    """
    Downloads large files that SimCADO needs to simulate MICADO
    """

    save_dir = os.path.join(__pkg_dir__, "data")
    fname = os.path.join(save_dir, "extras.dat")

    # check_replace = 0  ## unused (OC)
    if os.path.exists(fname):
        old_extras = ioascii.read(fname)
        # check_replace = 1   ## unused (OC)
    else:
        old_extras = ioascii.read("""
        filename                version         size    group
        PSF_POPPY.fits          20151103a       48MB    typical
        """)

    url = "http://www.univie.ac.at/simcado/data_ext/"
    new_extras = ioascii.read(download_file(url + "extras.dat"))

    for name, vers, size, group in new_extras:
        check_download = 1

        # does the file exist on the users disk?
        fname = os.path.join(__pkg_dir__, "data", name)
        if os.path.exists(fname):

            # is the new name in the old list of filenames
            if name in old_extras["filename"]:
                iname = np.where(old_extras["filename"] == name)[0][0]
                # print(iname, old_extras["version"][iname] == vers)

                # Are the versions the same?
                if vers == old_extras["version"][iname]:
                    check_download = 0

        if check_download:
            print("Downloading: " + name + "  Version: " + vers +
                  "  Size: " + size)
            download_file(url + name)
        else:
            print(name + " is already the latest version: " + vers)

    print("Finished downloading data for SimCADO")

    
def add_SED_to_simcado(file_in, file_out=None, lam_units="um"):
    """
    Adds the SED given in ``file_in`` to the SimCADO data directory
    
    Parameters
    ----------
    file_in : str
        path to the SED file. Can be either FITS or ASCII format with 2 columns
        Column 1 is the wavelength, column 2 is the flux
    file_out : str, optional
        Default is None. The file path to save the ASCII file. If ``None``, the SED 
        is saved to the SimCADO data directory i.e. to ``<utils.__pkg_dir__>/data/``
    lam_units : str, astropy.Units
        Units for the wavelength column, either as a string or as astropy units
        Default is [um]
    
    """
    
    file_name, file_ext = os.path.basename(file_in).split(".")
    
    if file_out is None:
        if "SED_" not in file_name:
            file_out = __pkg_dir__+"/data/SED_"+file_name+".dat"
        else: file_out = __pkg_dir__+"/data/"+file_name+".dat"
            
    if file_ext.lower() in "fits":
        data = fits.getdata(file_in)
        lam, val = data[data.columns[0].name], data[data.columns[1].name]
    else:
        lam, val = ioascii.read(file_in)[:2]

    lam = (lam * u.Unit(lam_units)).to(u.um)
    mask = (lam > 0.3*u.um) * (lam < 5.0*u.um) 

    np.savetxt(file_out, np.array((lam[mask], val[mask]), dtype=np.float32).T, 
               header="wavelength    value \n [um]         [flux]")

    
def airmass_to_zenith_dist(airmass):
    """
    returns zenith distance in degrees
    
    Z = arccos(1/X)
    """
    return np.rad2deg(np.arccos(1. / airmass))
    

def zentih_dist_to_airmass(zenith_dist):
    """
    ``zenith_dist`` is in degrees
    
    X = sec(Z)
    """
    return 1. / np.cos(np.deg2rad(zenith_dist))
   