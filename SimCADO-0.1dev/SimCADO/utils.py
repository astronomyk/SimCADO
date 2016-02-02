###############################################################################
# AtmosphereModel
#
# DESCRIPTION
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
#  read_config(config_file)
#  update_config(config_file, config_dict)
#  unify(x, unit, length=1)
#  parallactic_angle(ha, de, lat=-24.589167)
#  parallactic_angle_2(ha, de, lat=-24.589167)
#  moffat(r, alpha, beta)
#
#
#
#



import numpy as np
from astropy import units as u

## These functions are exported to the package
__all__ = ["read_config", "update_config", "unify", "parallactic_angle"]

## CHECK: Turn config into a class? (subclass of dict, if anything)
def read_config(config_file):
    '''Read SimCADO configuration file

    Input:
    -----
    A configuration file in the SExtractor format:
       'PARAMETER    Value    # Comment'

    Returns:
    -------
    A dictionary with keys 'PARAMETER' and values 'Value'.
    Note that the values of the dictionary are strings and will have 
    to be converted to the appropriate data type as they are needed. 
    '''

    import re
    
    # Read the file into a list of strings
    config_dict = dict()

    # remove lines that are all spaces or spaces + '#'
    # these are the regular expressions
    isempty = re.compile('^\s*$')
    iscomment = re.compile('^\s*#')

    with open(config_file, 'r') as fp:
        for line in fp:
            if isempty.match(line): continue
            if iscomment.match(line): continue
            line = line.rstrip()   # remove trailing \n
            try:
                (content, comment) = line.split('#', 1)
            except:
                content = line
                comment = ""
            param, value = content.split(None, 1)
            config_dict[param] = value.strip()

    return config_dict


def update_config(config_file, config_dict):
    '''Update SimCADO configuration

    Input:
    -----
    - A configuration file in the SExtractor format:
         'PARAMETER    Value    # Comment'
    - an existing configuration dictionary. 

    Returns:
    -------
    A dictionary with keys 'PARAMETER' and values 'Value'.
    Note that the values of the dictionary are strings and will have 
    to be converted to the appropriate data type as they are needed. 
    '''
    config_dict.update(read_config(config_file))

    return config_dict

def unify(x, unit, length=1):
    """Convert all types of input to an astropy array/unit pair"""

    print type(x)
    
    if type(x) == u.quantity.Quantity:
        if type(x.value) == np.ndarray:
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
    """Compute parallactic angle

    Parameters:
    ==========

    - ha, de: hour angle and declination of observation point
    - lat: latitude of observatory
    """

    # Convert observation point, pole and zenith to cartesian coordinates
    x_pole = [0, 0, -1]
    print x_pole
    
    poledist = (90. - lat)/180 * np.pi  # from north pole
    x_zenith = [0, np.sin(poledist), np.cos(poledist)]
    print x_zenith
    
    obsdist = (90. - de)/180 * np.pi
    ha = ha/180. * np.pi
    x_obs = [np.sin(obsdist) * np.sin(ha), np.sin(obsdist) * np.cos(ha),
             np.cos(obsdist)]
    print x_obs
    
    # normals to the great circles
    N_pole = np.cross(x_obs, x_pole)
    N_pole = N_pole/np.sqrt((N_pole**2).sum())
    N_zenith = np.cross(x_obs, x_zenith)
    N_zenith = N_zenith/np.sqrt((N_zenith**2).sum())

    # Angle between the great circles
    ## CHECK: exact definition, modulo 180 deg?
    parang = np.arccos((N_pole * N_zenith).sum()) * 180. / np.pi

    return parang

# From Filippenko1982
def parallactic_angle_2(ha, de, lat=-24.589167):
    hadeg = ha
    ha = ha/180. * np.pi
    de = de/180. * np.pi
    lat = lat/180. * np.pi

    sineta = np.sin(ha)*np.cos(lat)/np.sqrt(1. - (np.sin(lat) * np.sin(de) + np.cos(lat) * np.cos(de) * np.cos(ha))**2)

    eta =  np.arcsin(sineta) * 180./np.pi
    if hadeg >=0 and hadeg <=45:
        eta = 180. - eta
    elif hadeg < 0 and hadeg >=-45:
        eta = 180. + eta
    elif hadeg < -45:
        eta = - eta
    return eta

def moffat(r, alpha, beta):
    return (beta - 1)/(np.pi * alpha**2) * (1 + (r/alpha)**2)**(-beta)

	
	
	
