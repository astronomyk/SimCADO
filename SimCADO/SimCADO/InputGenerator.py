"""InputGenerator"""
################################################
# InputGenerator
#
# Put various functions here that generate input data for a simulation
#
# - FITS cubes with PSFs
# - FITS cubes with source data
# - ASCII lists with source data
# - Broadband filter curves, unity curves
# - Detector noise image/cubes
# - Dead pixel image
#

import os
import warnings
import numpy as np
from astropy.io import fits, ascii
import astropy.units as u

try:
    import SimCADO.utils as utils
    import SimCADO.LightObject as lo
except:
    import utils
    import LightObject as lo
    

def source_star(spec_type, distance=10*u.pc, filename=None, **kwargs):
    """
    create a Source object for a main sequence star. If filename is None, return
    the object. Otherwise save it to disk
    """
    # Must have the following properties:
    # ===================================
    # - lam       : LAM_MIN, LAM_MAX [, CDELT3, CRPIX3, CRVAL3, NAXIS3] 
    # - spectra   :
    # - x         : X_COL
    # - y         : Y_COL
    # - spec_ref  : REF_COL
    # - weight    : W_COL
    # - units*    : BUNIT
    # - pix_res*  : PIX_RES [, CDELT1]
    # - exptime*  : EXPTIME
    # - area*     : AREA

    params = {  "x"     :0
                "y"     :0
                "units" :"ph/(s m2)"
                "pix_res":0.004
                "exptime":None
                "area"  :None   }
    params.update(kwargs)
    
    lam, spectra = get_MS_spectra(spec_type, distance)
    x, y = [0], [0]
    spec_ref = [1]
    weight = [0]
    units = 
    pix_res = 
    

def source_from_mags():
    """
    Generate a LightObject FITS file from a list of magnitudes
    """
    pass




def poppy_eelt_psf_cube(lam_bin_centers, filename=None, **kwargs):
    """
    Generate a FITS file with E-ELT PSFs for a range of wavelengths by using the
    POPPY module - https://pythonhosted.org/poppy/
    
    Parameters:
    -----------
    lam_bin_centers: [um] the centres of each wavelength bin
    filename: [None] path name to where the FITS file should be saved. If not 
              given, a HDUList is returned
    
    Optional Parameters:
    --------------------
    pix_res: [arcsec] the angular resolution of the pixels. Default is 1 mas
    diameter_out: [m]
    diameter_in: [m]
    flatotflat: [m]
    gap: [m]
    n_spiders: [m]
    size: [int]
    oversample: [int]
    clobber: [True/False]
    """
    
    try:
        import poppy
    except:
        raise ValueError("Please install poppy \n >>sudo pip3 install poppy")
        return
    
    params = {  "diameter_out"  :37,
                "diameter_in"   :5.5,
                "flattoflat"    :1.45,
                "gap"           :0.004,
                "n_spiders"     :6,
                "pix_res"       :0.001,
                "size"          :255,
                "oversample"    :1,
                "clobber"       :True   }
    params.update(kwargs)
    
    rings = int(0.65 * params["diameter_out"] / params["flattoflat"])
    m1 = poppy.MultiHexagonAperture(rings=rings, 
                                    flattoflat=params["flattoflat"], 
                                    gap=params["gap"])
    pri = poppy.CircularAperture(radius=params["diameter_out"]/2)
    sec = poppy.SecondaryObscuration(secondary_radius=params["diameter_out"]/2,
                                     n_supports=params["n_spiders"], 
                                     support_width=0.5)
    eelt = poppy.CompoundAnalyticOptic( opticslist=[m1, pri, sec], name='E-ELT')

    osys = poppy.OpticalSystem()
    osys.addPupil(eelt)
    osys.addDetector(pixelscale=params["pix_res"],
                     fov_arcsec=params["pix_res"] * params["size"], 
                     oversample=params["oversample"])

    psfHDU = [osys.calcPSF(lam * 1E-6)[0] for lam in lam_bin_centers]
    hdulist = fits.HDUList(psfHDU)
    
    if filename is None:
        return hdulist
    else: 
        hdulist.writeto(filename, clobber=params["clobber"])

        
        
    
################################################################################    
#              Stellar parameters from mass for InputGenerator                 #
################################################################################

# The following functions are to help determine various stellar parameters based
# on the mass of the star. The is so that a cluster of main sequence stars can 
# be generated according to masses in an IMF.
    
    
        
        
        
def mass_to_temp(mass):
    
    mass = mass.value if type(mass) == u.quantity.Quantity else mass
    
    f = np.array([ 0.02651303, -0.05307791, -0.10533279,  0.1843677 ,  0.5460582 , 3.74004826])
    logM = np.log10(mass)
    logT = np.polyval(f, logM)
    
    return 10**logT * u.K
        
    
def mass_to_flux5556A(mass, distance=10*u.pc):    
    """
    Calculate a rough approximation of the 5556 Angstrom flux of a star. 
    Units returned are ph/s/m2/um
    
    Parameters
    ==========
    - mass : float 
        [solar masses] the mass of the star
    - distance : float, optional 
        [parsec]the distance to the star
    
    F = F0 * 10**(-0.4 * f(log10(temp)))
    
    The difference in luminosity is based on the difference in absolute 
    magnitude, Mv, compared to Vega (@ 10pc Mv=0) and Mv is proportional to 
    log10(Temp), which we get from the mass.
    
    The value for F0 is taken to be 996 ph/s/cm2/A (or 99.6E9 ph/s/m2/um)
    for a Mv = 0 star
    from "The Observation and Analysis of Stellar Photospheres - D. Gray"
    
    F_5556A = f(mass)
    """
    
    mass = mass.value          if type(mass) == u.quantity.Quantity     else mass
    distance = distance * u.pc if type(distance) != u.quantity.Quantity else distance
    temp = mass_to_temp(mass)
    
    f = np.array([-10.57984516,  139.88537641, -624.14454692,  935.92676894])
    
    if type(temp) == u.quantity.Quantity:
        logT = np.log10(temp.value) 
    else: 
        np.log10(temp)
    Mv = np.polyval(f, logT)

    # janskys = (3580 * u.Jy * 10**(-0.4*Mv) * (10*u.pc / distance)**2).to(u.Jy)
    photons = (996 * u.Unit("ph/(s cm2 Angstrom)") * 10**(-0.4*Mv) * \
                                            (10*u.pc / distance)**2)
    return photons.to(u.ph/u.s/u.m**2/u.um).value
        
    
def mass_to_lifetime(mass):
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
    
    
def mass_to_spec_type(mass):
    """
    Returns the spectral type of a main sequence star with a given mass, 
    according to a 6th order polynomical best fit to the data set: 
        
        f = np.polyfit(mass, n, 6)
        
    where n is the spectral type in numerical form: O5 = 5, A0 = 20, G2 = 42,etc
    """
    d = {0:"o", 1:"b", 2:"a", 3:"f", 4:"g", 5:"k", 6:"m", 7:"l", 8:"t"}
    #f = np.polyfit(logM, logN, 6)
    f = np.array([ 0.02339193, -0.15822858,  0.0989074 ,  0.38510294, 
                  -0.30572525, -0.61230356,  1.62685421])
    
    logM = np.log10(mass)
    n = np.array(10**np.polyval(f, logM), dtype=int)
    spec_type = d[n//10] + str(n%10)
    return spec_type

    
def spec_type_to_mass(spec_type):
    """
    Takes a spectral type of a main sequence star and returns the mass, 
    according to a 8th order polynomical best fit to the data set: 
        
        f = np.polyfit(n, mass, 8)
        
    where n is the spectral type in numerical form: O5 = 5, A0 = 20, G2 = 42,etc
    """
    d = {"o":0, "b":1, "a":2, "f":3, "g":4, "k":5, "m":6, "l":7, "t":8}
    #f = np.polyfit(np.log10(n), np.log10(mass), 8)
    f = np.array([ -15.63725949,  114.5328907 , -346.14762635,  555.10413625,
                  -503.67053477,  253.88948443,  -65.56109921,    6.66706868,
                     1.98858127])
    
    if type(spec_type) == list:
        n = [float(d[SpT[0].lower()])*10 + float(SpT[1]) for SpT in spec_type]
    else:
        n = float(d[spec_type[0].lower()])*10 + float(spec_type[1])
    
    logN = np.log10(n)
    mass = np.round(10**np.polyval(f, logN), 3)   
    return mass    
    
    
def spec_type_to_n(spec_type):    
    """
    Convert a spectral type to numerical type: O5 = 5, A0V = 20, G2 = 42, etc
    O:0, B:1, A:2, F:3, G:4, K:5, M:6, L:7, T:8
    """
    d = {"o":0, "b":1, "a":2, "f":3, "g":4, "k":5, "m":6, "l":7, "t":8}
    n = float(d[spec_type[0].lower()])*10 + float(spec_type[1])
    return int(n)


def n_to_spec_type(n):    
    """
    Convert a numerical type to spectral type : O5 = 5, A0V = 20, G2 = 42, etc
    O:0, B:1, A:2, F:3, G:4, K:5, M:6, L:7, T:8
    """
    d = {0:"o", 1:"b", 2:"a", 3:"f", 4:"g", 5:"k", 6:"m", 7:"l", 8:"t"}
    spec_type = d[n//10] + str(n%10) + "v"
    return spec_type    
    
    
def mass_to_n_class(mass):
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
    
    
def mass_to_pickles(mass, pickles_dir):
    """
    Not all main sequence spectral types are in the Pickles library. This 
    function looks at what is available and then assigns a spectrum to each star
    based on its mass
    """
    n = n_class_from_mass(mass)
    stars   = [fname[2:4] for fname in os.listdir(pickles_dir) if "uk" in fname and "v" in fname[4]]
    avail   = [int(str(["o","b","a","f","g","k","m"].index(f[0]))+f[1]) for f in stars]
    dn      = [(avail-i)[abs(avail-i).argmin()] for i in n]

    x,y = (n+dn)/10, (n+dn)%10
    pickle_types = [["o","b","a","f","g","k","m"][x[i]]+str(y[i])+"v" for i in range(len(n))]
    
    return pickle_types
    
    
def nearest_pickle_type(spec_type):
    """
    Get the closest spectrum out of the incomplete catalogue from Pickles 1998
    http://adsabs.harvard.edu/abs/1998PASP..110..863P
    """
    keys = ['o5v', 'o9v', 
            'b0v', 'b1v', 'b3v', 'b8v', 'b9v', 
            'a0v', 'a2v', 'a3v', 'a5v', 'a7v',
            'f0v', 'f2v', 'f5v', 'f6v', 'f8v', 
            'g0v', 'g2v', 'g5v', 'g8v', 
            'k0v', 'k2v', 'k3v', 'k4v', 'k5v', 'k7v', 
            'm0v', 'm1v', 'm2v', 'm3v', 'm4v', 'm5v', 'm6v']
    n = np.array([spec_type_to_n(i) for i in keys if i[1].isdigit()])
    
    if type(spec_type) == str:
        if spec_type[1].isdigit():
            m = spec_type_to_n(spec_type)
        else:
            raise ValueError(spec_type+" isn't a main sequence star")
    else:
        m = int(spec_type)

    pick_type = n_to_spec_type(n[np.argmin(np.abs(n-m))])
    return pick_type
    
  
def get_MS_spectra(spec_type, distance=10*u.pc, pickles_table=None):
    """
    Pull in the emission curve for a specific star
    
    Parameters
    ==========
    - spec_type : str, [str]
        the spectral class(es) of the main sequence star(s), e.g. A0, G2V
    - distance : float, [float], optional
        [parsec] the distance(s) to said star
    - star_catalogue : str, optional
        path name to an ASCII table with spectra of MS stars normalised to
        lam = 5556 Angstrom. Default is "../data/EC_pickles_MS.dat".
        
    Notes:
        It is assumed that the wavelength column of any pickles_table is in [um]
        Units of the returned spectra are [ph/s/m2]
    """   
    if pickles_table is None:
        pickles_table = ascii.read("../data/EC_pickles_MS.dat")
    cat = pickles_table  
    
    if type(spec_type) == list:
        pick_type = [nearest_pickle_type(SpT) for SpT in spec_type]
    else:
        pick_type = nearest_pickle_type(spec_type)
        
    # Flux is in units of [ph/s/m2/um]. dlam is in units of [um]    
    mass = spec_type_to_mass(pick_type)
    flux = mass_to_flux5556A(mass, distance)
    
    lam = cat[cat.colnames[0]].data
    dlam = np.append(lam[1:] - lam[:-1], lam[-1] - lam[-2])
    
      
    if type(spec_type) == list:
        spectra = []
        for i in range(len(pick_type)):
            tmp = cat[pick_type[i][:2].lower()+"v"].data
            spectra += [tmp * dlam * flux[i]]
    else:
        spectra = cat[spec_type[:2].lower()+"v"].data
        spectra = spectra * dlam * flux
    
    return lam, spectra 

    