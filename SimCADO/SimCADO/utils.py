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

## These functions are exported to the package
__all__ = ["read_config", "update_config", "unify", "parallactic_angle", "atmospheric_refraction"]

## CHECK: Turn config into a class? (subclass of dict, if anything)
def read_config(config_file):
	'''Read SimCADO configuration file

	Input:
	-----
	A configuration file in the SExtractor format:
	   'PARAMETER	 Value	  # Comment'

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
			try: config_dict[param] = float(value.strip())
			except: config_dict[param] = value.strip()

	return config_dict


def update_config(config_file, config_dict):
	'''Update SimCADO configuration

	Input:
	-----
	- A configuration file in the SExtractor format:
		 'PARAMETER	   Value	# Comment'
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

	print(type(x))
	
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
	print(x_pole)
	
	poledist = (90. - lat)/180 * np.pi	# from north pole
	x_zenith = [0, np.sin(poledist), np.cos(poledist)]
	print(x_zenith)
	
	obsdist = (90. - de)/180 * np.pi
	ha = ha/180. * np.pi
	x_obs = [np.sin(obsdist) * np.sin(ha), np.sin(obsdist) * np.cos(ha),
			 np.cos(obsdist)]
	print(x_obs)
	
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


def atmospheric_refraction(lam, z0=60, temp=0, rel_hum=60, pres=750, lat=-24.5, h=3064):
	"""
	Work out the angle [arcsec] difference between where the star is above and below the atmosphere
	[See Stone 1996 and the review by S. Pedraz - http://www.caha.es/newsletter/news03b/pedraz/newslet.html]
	# Input
	# lam = [micron] wavelength 
	# z0 = [degrees] zenith distance 
	# temp = [deg C] ground temperature 
	# rel_hum = [%] relative humidity 
	# pres = [millibar] air pressure 
	# 
	# Constants (which make almost no difference)
	# lat = [degrees] latitude  - set for Cerro Armazones: 24.5 deg South
	# h = [meters] height above sea level  - 3064 m
	#
	# Output
	# ang = [arcsec] angle between real position and refracted position 
	"""
	#### need T, P, RH for Ps, Pw Pa
	T = 273.15 + temp

	Ps = -10474. + 116.43*T - 0.43284*T**2 + 0.0005384*T**3
	Pw = rel_hum/100.*Ps
	Pa = pres - Pw

	#### need n0 for gamma
	sig = 1./lam
	Da = (Pa/T) * (1. + Pa*(57.9E-8 - 0.0009325/T + 0.25844/T**2))
	Dw = (Pw/T) * (1. + Pw*(1. + 3.7E-4*Pw)*(-2.37321E-3 + 2.23366/T - 710.792/T**2 + 77514.1/T**3))

	na = Da*(2371.34 + 683939.7/(130.-sig**2) + 4547.3/(38.9-sig**2))
	nw = Dw*(6487.31 + 58.058*sig**2 - 0.7115*sig**4 + 0.08851*sig**6)
	n0 = 1E-8*(na + nw) + 1.

	#### need gamma, kappa and beta for R
	g = n0 - 1.
	b = 0.001254 * (273.15 + temp) / 273.15
	#k = 1
	k = 1. + 0.005302*np.sin(lat/57.29578)**2 - 0.00000583*np.sin(2.*lat/57.29578)**2 - 0.000000315*h
	
	R = k*g*(1 - b) * np.tan(z0/57.29578) - k*g*(b - g/2.) * np.tan(z0/57.29578)**3

	# the refraction is the side of a triangle, although the triangle side makes an arc across the sky. 
	# We want the angle that this triangle side subtends
	# Using the small angle approximation (which is in radians), we can get the angle of refraction
	
	ang = R * 3600 * 57.29578
	
	return ang	

    
def poissonify(self, arr):
    """ 
    Add a realisation of the poisson process to the array 'arr'. 

    Keywords:
    - arr: 
    """
    return np.random.poisson(arr).astype(np.float32)
    
    
    
    
################################################################################    
#              Stellar parameters from mass for InputGenerator                 #
################################################################################

# The following functions are to help determine various stellar parameters based
# on the mass of the star. The is so that a cluster of main sequence stars can 
# be generated according to masses in an IMF.
    
    
        
        
        
def temp_from_mass(mass):
    
    mass = mass.value if type(mass) == u.quantity.Quantity else mass
    
    f = np.array([ 0.02651303, -0.05307791, -0.10533279,  0.1843677 ,  0.5460582 , 3.74004826])
    logM = np.log10(mass)
    logT = np.polyval(f, logM)
    
    return 10**logT * u.K
        
    
def flux5556A_from_mass(mass, distance=10*u.pc):    
	"""
	F = F0 * 10**(-0.4 * f(log10(temp)))
	
	The difference in luminosity is based on the difference in absolute 
    magnitude, Mv, compared to Vega (@ 10pc Mv=0) and Mv is proportional to 
    log10(Temp), which we get from the mass.
	
	The value for F0 is taken to be 3580 Jy for a Mv = 0 star
	This is based on the value for Vega, at 7.7pc with Mv(vega) = 0.58
	Therefore if Vega were at 10pc, it would have a Mv = 0.
    
    F_5556A = f(mass)
	"""
	
	mass = mass.value if type(mass) == u.quantity.Quantity else mass
	temp = temp_from_mass(mass)
	
	f = np.array([-10.57984516,  139.88537641, -624.14454692,  935.92676894])
	
	if type(temp) == u.quantity.Quantity:
        logT = np.log10(temp.value) 
    else: 
        np.log10(temp)
	Mv = np.polyval(f, logT)
	
	return (3580 * u.Jy * 10**(-0.4*Mv) * (10*u.pc / distance)**2).to(u.Jy)
    
    
def lifetime_from_mass(mass):
    """
    Calculate the lifetime of a main sequence star based on a best fit
    approximation to a plot of lifetimes vs masses.
    
    t = f(mass)
    """
    mass = mass.value   if type(mass) == u.quantity.Quantity    else mass
    
    f = np.array([ 0.24897294, -0.09842223, -2.47114954,  9.89215499])
    
    logM = np.log10(mass)
    logt = np.polyval(f, logM)
    
    return 10**logt * u.yr
    
    
def n_class_from_mass(mass):
	"""
    I have given spectral types a number based (O=0, A=2, M=6, etc) so that they
    can be quantified. The following determines the numerical main sequence 
    spectral type based on the stars mass.
    
    n(spec_type) = f(mass)
    """
	mass = mass.value if type(mass) == u.quantity.Quantity else mass
	
	f = np.array([-0.14376899,  0.13219846,  0.37555566, -0.31116164, -0.59572498, 1.62977733])
	logM = np.log10(mass)
	logN = np.polyval(f, logM)
	
	n = np.asarray(10**logN, dtype=int)
	return n
	
	
def pickles_from_n(mass, pickles_dir):
	"""
    Not all main sequence spectral types are in the Pickles library. This 
    function looks at what is available and then assigns a spectrum to each star
    based on its mass
    """
    n = n_class_from_mass(mass)
	stars 	= [fname[2:4] for fname in os.listdir(pickles_dir) if "uk" in fname and "v" in fname[4]]
	avail	= [int(str(["o","b","a","f","g","k","m"].index(f[0]))+f[1]) for f in stars]
	dn      = [(avail-i)[abs(avail-i).argmin()] for i in n]

	x,y = (n+dn)/10, (n+dn)%10
	pickle_types = [["o","b","a","f","g","k","m"][x[i]]+str(y[i])+"v" for i in range(len(n))]
	
	return pickle_types
    