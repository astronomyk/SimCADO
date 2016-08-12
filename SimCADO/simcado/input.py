"""opticsGenerator"""
################################################
# opticsGenerator
#
# Put various functions here that generate optics data for a simulation
#
# - FITS cubes with PSFs
# - FITS cubes with source data
# - ASCII lists with source data
# - Broadband filter curves, unity curves
# - detector noise image/cubes
# - Dead pixel image
#

from datetime import datetime as dt
import numpy as np
from astropy.io import fits
from astropy.io import ascii as ioascii  # 'ascii' redefines built-in
import astropy.units as u

try:
    import simcado.source as lo
    import simcado.spectral as sc
    import simcado.psf as psf
except ImportError:
    import source as lo
    import spectral as sc
    import psf as psf


__all__ = ["diff_limited_eelt_psf_cube", "poppy_eelt_psf_cube",
           "grid_of_stars", "stellar_emission_curve", "source_1E4_Msun_cluster",
           "source_from_stars"]


################################################################################
#                       Generate Instrument optics Data                         #
################################################################################

def diff_limited_eelt_psf_cube(lam_bin_centers, filename=None, **kwargs):
    """
    Generate a FITS file with a diffraction limited E-ELT PSFs for a range of
    wavelengths by using psf.AiryPSFCube()

    Parameters:
    -----------
    lam_bin_centers: [um] the centres of each wavelength bin
    filename: [None] path name to where the FITS file should be saved. If not
              given, a psf object is returned

    Optional Parameters:
    --------------------
    pix_res: [arcsec] the angular resolution of the pixels. Default is 1 mas
    diameter: [m]
    size: [int]
    clobber: [True/False]
    """

    params = {"diameter"      :37,
              "pix_res"       :0.001,
              "size"          :255,
              "clobber"       :True}
    params.update(kwargs)

    my_psf = psf.AiryPSFCube(lam_bin_centers, **params)

    if filename is None:
        return my_psf
    else:
        my_psf.export_to_fits(filename, clobber=params["clobber"])


def poppy_eelt_psf_cube(lam_bin_centers, filename=None, **kwargs):
    """
    Generate a FITS file with E-ELT PSFs for a range of wavelengths by using the
    POPPY module - https://pythonhosted.org/poppy/

    Parameters:
    -----------
    lam_bin_centers: [um] the centres of each wavelength bin in micron
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
    except ImportError:
        print("Please install poppy \n >>sudo pip3 install poppy")
        return

    params = {"diameter_out"  :37,
              "diameter_in"   :5.5,
              "flattoflat"    :1.45,
              "gap"           :0.004,
              "n_spiders"     :6,
              "pix_res"       :0.001,
              "size"          :255,
              "oversample"    :1,
              "clobber"       :True,
              "cpus"          :-1}
    params.update(kwargs)

    # Create the Optical Train
    rings = int(0.65 * params["diameter_out"] / params["flattoflat"])

    m1 = poppy.MultiHexagonAperture(rings=rings,
                                    flattoflat=params["flattoflat"],
                                    gap=params["gap"])
    pri = poppy.CircularAperture(radius=params["diameter_out"]/2)
    sec = poppy.SecondaryObscuration(secondary_radius=params["diameter_in"]/2,
                                     n_supports=params["n_spiders"],
                                     support_width=0.5)
    eelt = poppy.CompoundAnalyticOptic(opticslist=[m1, pri, sec], name='E-ELT')

    # Create a "detector"
    osys = poppy.OpticalSystem()
    osys.addPupil(eelt)
    osys.add_detector(pixelscale=params["pix_res"],
                      fov_arcsec=params["pix_res"] * params["size"],
                      oversample=params["oversample"])


    # list of wavelengths for which I should generate PSFs
    lam_bin_centers = np.array(lam_bin_centers)

    # Generate the PSFs in multiple threads
    try:
        import multiprocessing as mp
        if params["cpus"] < 0:
            num_cpu = max(1, mp.cpu_count() - 1)
        else:
            num_cpu = min(params["cpus"], mp.cpu_count() - 1)

        pool = mp.Pool(processes=num_cpu)
        pond = [pool.apply_async(_get_poppy_psf, (osys, lam)) \
                                                    for lam in lam_bin_centers]
        psf_hdu = [duck.get() for duck in pond]

    except (ImportError, mp.ProcessError): # any other exception to catch? (OC)
        print("Not possible to use multiple threads")
        psf_hdu = []
        for lam in lam_bin_centers:
            psf_hdu += [_get_poppy_psf(osys, lam)]
            print(lam)

    for psfi in psf_hdu:  # was psf, redefined module! (OC)
        psfi.data = psfi.data.astype(np.float32)
        psfi.header["CDELT1"] = (params["pix_res"],
                                 "[arcsec] - Pixel resolution")
        psfi.header["CDELT2"] = (params["pix_res"],
                                 "[arcsec] - Pixel resolution")
        psfi.header["PSF_TYPE"] = ("POPPY", "Generated by Poppy")
        # Convert wavelength in header to [um] - originally in [m]
        psfi.header["WAVE0"] = (psfi.header["WAVE0"]*1E6,
                                "Wavelength of PSF [um]")

    if len(psf_hdu) > 1:
        psf_hdu = [psf_hdu[0]] + \
                 [fits.ImageHDU(i.data, i.header) for i in psf_hdu[1:]]

    hdulist = fits.HDUList(psf_hdu)

    if filename is None:
        return hdulist
    else:
        hdulist.writeto(filename, clobber=params["clobber"])


def _get_poppy_psf(osys, lam):
    """
    A self contained function for "reading out" the detector, so that we can
    use multi-threading, if available (i.e. if it isn't a windows machine)
    """
    t = dt.now()
    print(lam, str(t.hour)+":"+str(t.minute)+":"+str(t.second))
    return osys.calcPSF(lam * 1E-6)[0]




################################################################################
#                       Generate some Source objects                           #
################################################################################

def grid_of_stars(mag_range, spacing=0.5, spec_type="a0v", filename=None):
    """

    spacing : float, optional
        [arcsec] distance between each star. Default = 0.5 arcsec
    """


    n = len(mag_range)
    side_len = int(np.sqrt(n)) + (np.sqrt(n) % 1 > 0)
    print(side_len)
    x = spacing * (np.arange(n) % side_len - (side_len - 1) / 2)
    y = spacing * (np.arange(n)// side_len - (side_len - 1) / 2)

    return source_from_stars(spec_type.lower()*len(mag_range),
                             x, y,
                             dist_mod=mag_range,
                             filename=filename)


def stellar_emission_curve(spec_type, distance=None, dist_mod=None,
                           in_photons=True, **kwargs):
    """
    Get an emission curve for a certain type of star
    !! TODO scale the emission curve to the right magnitude
    """
    if distance is None:
        if dist_mod is None:
            distance = 10*u.pc
        else:
            distance = 10**(1. + 0.2 * dist_mod)

    lam, spec = get_MS_spectra(spec_type, distance=distance,
                               in_photons=in_photons)
    lam_res = np.median(lam[1:]-lam[:-1])

    #lam_central = {"u":0.36, "b":0.44, "v":0.55, "r":0.64, "i":0.79,
    #               "j":1.26, "h":1.60, "k":2.22, "ks":2.16}
    # not done yet
    #a0v_rel_flux = {"u":0.36, "b":0.44, "v":0.55, "r":0.64, "i":0.79,
    #                "j":1.26, "h":1.60, "k":2.22, "ks":2.16}

    return sc.EmissionCurve(lam=lam, val=spec, lam_res=lam_res,
                            units="ph/(s m2)", **kwargs)


def source_1E4_Msun_cluster(distance=40000, half_light_radius=1, filename=None,
                            **kwargs):
    """
    Generate a source object for a 10^4 solar mass cluster

    Parameters
    ==========
    distance : float
        [pc] distance to the cluster
    half_light_radius : float
        [pc] half light radius of the cluster
    filename : str, optional
        Path to where the Source object FITS file is saved.  If filename is None
        the object is returned.

    **kwargs
        {units  : "ph/(s m2)", pix_res : 0.004, exptime : None, area : None}

    Returns a Source object
    """
    params = {"pix_res" :0.004}
    params.update(kwargs)

    masses = np.loadtxt("../data/IMF_1E4.dat")
    spec_type = [mass_to_spec_type(mass) for mass in masses]

    distance *= u.pc
    half_light_radius *= u.pc
    hwhm = (half_light_radius/distance*u.rad).to(u.arcsec).value
    sig = hwhm / 1.175

    x = np.random.normal(0, sig, len(masses))
    y = np.random.normal(0, sig, len(masses))

    return source_from_stars(spec_type, x, y, distance, filename=filename,
                             **kwargs)



def source_from_stars(spec_type, x, y, distance=None, dist_mod=None,
                      filename=None, **kwargs):
    """
    Create a source.Source object for a main sequence star.

    Parameters
    ==========
    spec_type : str
        Any Main sequence spectral type
    x : float, array-like
        [arcsec]
    y : float, array-like
        [arcsec]
    distance : float, array-like, optional
        Distance to the stars. Default is 10pc
    dist_mod : float, array-like, optional
        Distance modulus to the stars. If distance is given, dist_mod is ignored
    filename : str, optional
        Path to where the Source object FITS file is saved.  If filename is None
        the object is return.

    **kwargs
        {units  : "ph/(s m2)", pix_res : 0.004, exptime : None, area : None}
    """
    # Must have the following properties:
    # ===================================
    # - lam       : LAM_MIN, LAM_MAX [, CDELT3, CRPIX3, CRVAL3, NAXIS3]
    # - spectra   :
    # - x         : X_COL
    # - y         : Y_COL
    # - spec_ref  : REF_COL
    # - weight    : W_COL

    # **kwargs
    # - units*    : BUNIT
    # - pix_res*  : PIX_RES [, CDELT1]
    # - exptime*  : EXPTIME
    # - area*     : AREA


    # The default **kwargs are based off the units. Default it is ph/(s m2),
    # so the area and exptime kwargs are set to None.
    # We are looking at point sources here so pix_res plays no role. It is only
    # needed for the Source object
    params = {"units"  : "ph/(s m2)",
              "pix_res": 0.004,
              "exptime": None,
              "area"   : None}
    params.update(kwargs)
    is_list = isinstance(spec_type, (list, tuple, np.ndarray))

    # Work out which spectral types are unique
    spec_type = [spec_type] if not is_list else spec_type
    unique_type = list(np.unique(spec_type))

    # Pull in the spectra of the unique types, all for a distance of 10pc
    ref = [unique_type.index(i) for i in spec_type]
    lam, spectra = get_MS_spectra(unique_type)
    x = [x] if not is_list else x
    y = [y] if not is_list else y

    if distance is None:
        if dist_mod is None:
            distance = 10.
        else:
            if isinstance(dist_mod, list):
                dist_mod = np.array(dist_mod)
            distance = 10**(1. + 0.2 * dist_mod)
    elif isinstance(distance, u.Quantity):
        distance = distance.to(u.pc).value
    elif isinstance(distance, list):
        distance = np.asarray(distance)

    # Weight the spectra according to distance
    weight = (10. / distance)**2
    if not isinstance(weight, np.ndarray):
        weight = np.asarray([weight]*len(spec_type))

    # Create a source.Source object and write it to a file
    obj = lo.Source(lam=lam, spectra=spectra, x=x, y=y, ref=ref, weight=weight,
                    **params)
    if filename is not None:
        obj.write(filename)
        print("Source object written to "+filename)
    else:
        return obj


def source_from_mags():
    """
    Generate a source FITS file from a list of magnitudes
    """
    pass


def get_MS_spectra(spec_type, distance=10*u.pc, pickles_table=None, raw=False,
                   in_photons=True):
    """
    Pull in the emission curve for a specific star

    Parameters
    ==========
    - spec_type : str, [str]
        the spectral class(es) of the main sequence star(s), e.g. A0, G2V
    - distance : float, [float], optional
        [parsec] the distance(s) to the star(s)
    - star_catalogue : str, optional
        path name to an ASCII table with spectra of MS stars normalised to
        lam = 5556 Angstrom. Default is "../data/EC_pickles_MS.dat".
    - raw : bool, optional
        if True, only the raw spectrum (normalised to 1 at lam=5556A) is returned
    - in_photons : bool, optional
        returns the spectrum in [ph/s/m2/um] if True. If False, []
    Notes:
        It is assumed that the wavelength column of any pickles_table is in [um]
        Units of the returned spectra are [ph/s/m2]
    """
    if pickles_table is None:
        pickles_table = ioascii.read("../data/EC_pickles_MS.dat")
    cat = pickles_table

    if isinstance(spec_type, (list, tuple)):
        pick_type = [nearest_pickle_type(SpT) for SpT in spec_type]
    else:
        pick_type = nearest_pickle_type(spec_type)

    lam = np.round(cat[cat.colnames[0]].data, 7)
    dlam = np.append(lam[1:] - lam[:-1], lam[-1] - lam[-2])

    # flux is in units of [ph/s/m2/um]. dlam is in units of [um]
    # flux is adjusted for distance
    mass = spec_type_to_mass(pick_type)
    flux = mass_to_flux5556A(mass, distance, in_photons=in_photons)

    if isinstance(spec_type, (list, tuple)):
        spectra = []
        for i, thispick in enumerate(pick_type):
            tmp = cat[thispick[:2].lower()+"v"].data
            if not raw:
                spectra += [tmp * dlam * flux[i]]
            else:
                spectra += [tmp]
    else:
        spectra = cat[pick_type[:2].lower()+"v"].data
        if not raw:
            spectra = spectra * dlam * flux

    if in_photons:
        spectra /= 0.5556*u.um / lam

    # units of the spectra being returned are [ph/s/m2] for each wavelength bin,
    # adjusted for distance
    return np.asarray(lam), np.asarray(spectra)



################################################################################
#              Stellar parameters from mass for opticsGenerator                 #
################################################################################

# The following functions are to help determine various stellar parameters based
# on the mass of the star. The is so that a cluster of main sequence stars can
# be generated according to masses in an IMF.


def mass_to_Mv(mass):
    """
    Calculate the V-band absolute magnitude of a star based on mass
    """
    ## provide reference for this
    mass = mass.value if isinstance(mass, u.quantity.Quantity) else mass

    f = np.array([ -0.40960919,   1.28000053,  -0.43203104,  -1.48032921,
                    2.74563706,  -4.61182633,  -0.38868252,   8.58781463,
                  -12.21570688,   4.90158585])
    Mv = np.polyval(f, np.log10(mass))

    return Mv



def mass_to_temp(mass):
    """
    Calculate an approximate surface temperature based on mass
    """
    ## provide reference for this
    mass = mass.value if isinstance(mass, u.quantity.Quantity) else mass

    f = np.array([ 0.02651303, -0.05307791, -0.10533279,  0.1843677 ,
                   0.5460582 , 3.74004826])
    logM = np.log10(mass)
    logT = np.polyval(f, logM)

    return 10**logT * u.K


def mass_to_flux5556A(mass, distance=10*u.pc, in_photons=True):
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
    magnitude, Mv, compared to Vega (@ 10pc Mv=0)
    Mv is proportional to log10(Temp), which we get from the mass.

    The value for F0 is taken to be 996 ph/s/cm2/A (or 99.6E9 ph/s/m2/um)
    for a Mv = 0 star at 10pc
    from "The Observation and Analysis of Stellar Photospheres - D. Gray"

    F_5556A = f(mass)
    """

    if isinstance(mass, u.quantity.Quantity):
        mass = mass.value
    if not isinstance(distance, u.quantity.Quantity):
        distance = distance * u.pc


    # Interestingly the round about way, i.e. fitting mass --> temp --> Mv
    # is more accurate than the direct fit - mass --> Mv
    temp = mass_to_temp(mass)

    f = np.array([-10.57984516, 139.88537641, -624.14454692, 935.92676894])
    logT = np.log10(temp.value)

    #Mv = np.polyval(f, logT)  # What does that do? (OC)
    Mv = mass_to_Mv(mass) - 1.1

    # janskys = (3580 * u.Jy * 10**(-0.4*Mv) * (10*u.pc / distance)**2).to(u.Jy)
    if in_photons:
        photons = (996 * u.Unit("ph/(s cm2 Angstrom)") * 10**(-0.4*Mv) * \
                                                        (10*u.pc / distance)**2)
        return photons.to(u.ph/u.s/u.m**2/u.um).value

    else:
        ergs = (363.1E-11 * u.Unit("erg/(s cm2 Angstrom)") * 10**(-0.4*Mv) * \
                                                        (10*u.pc / distance)**2)
        return ergs.to(u.erg/u.s/u.m**2/u.um).value




def mass_to_lifetime(mass):
    """
    Calculate the lifetime of a main sequence star based on a best fit
    approximation to a plot of lifetimes vs masses.

    t = f(mass)
    """
    if isinstance(mass, u.quantity.Quantity):
        mass = mass.value
    f = np.array([0.24897294, -0.09842223, -2.47114954, 9.89215499])

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
    ## Needs a reference (OC)
    d = {0:"o", 1:"b", 2:"a", 3:"f", 4:"g", 5:"k", 6:"m", 7:"l", 8:"t"}
    #f = np.polyfit(logM, logN, 6)
    f = np.array([ 0.02339193, -0.15822858,  0.0989074 ,  0.38510294,
                  -0.30572525, -0.61230356,  1.62685421])

    logM = np.log10(mass)
    n = np.array(10**np.polyval(f, logM), dtype=int)
    if n > 89:           # Err, what? n is supposed to be an array (OC)
        n = 89
    spec_type = d[n//10] + str(n%10)
    return spec_type


def mass_to_n_class(mass):
    """
    I have given spectral types a number based (O=0, A=2, M=6, etc) so that they
    can be quantified. The following determines the numerical main sequence
    spectral type based on the stars mass.

    n(spec_type) = f(mass)
    """
    if isinstance(mass, u.quantity.Quantity):
        mass = mass.value

    f = np.array([-0.14376899,  0.13219846,  0.37555566, -0.31116164,
                  -0.59572498, 1.62977733])
    logM = np.log10(mass)
    logN = np.polyval(f, logM)

    n = np.asarray(10**logN, dtype=int)
    return n


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

    if isinstance(spec_type, list):
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
    if not isinstance(spec_type, list):
        spec_type = [spec_type]
    d = {"o":0, "b":1, "a":2, "f":3, "g":4, "k":5, "m":6, "l":7, "t":8}
    n = [int(float(d[SpT[0].lower()])*10 + float(SpT[1])) for SpT in spec_type]
    return n


def n_to_spec_type(n):
    """
    Convert a numerical type to spectral type : O5 = 5, A0V = 20, G2 = 42, etc
    O:0, B:1, A:2, F:3, G:4, K:5, M:6, L:7, T:8
    """
    n = [n] if not hasattr(n, "__iter__") else n
    d = {0:"o", 1:"b", 2:"a", 3:"f", 4:"g", 5:"k", 6:"m", 7:"l", 8:"t"}
    spec_type = [d[i//10] + str(i%10) + "v" for i in n]
    return spec_type

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

    # n are the numerical spectral types, m is the numerical type of "spec_type"
    # find when n is closest to n, then convert that n to a string spectral type
    if not isinstance(spec_type, (list, tuple, np.array)):
        if isinstance(spec_type, (str, np.str_)):
            if spec_type[1].isdigit():
                m = spec_type_to_n(spec_type)
            else:
                raise ValueError(spec_type+" isn't a main sequence star")
        else:
            m = int(spec_type)

        pick_type = n_to_spec_type(n[np.argmin(np.abs(n-m))])
        return pick_type[0]
    else:
        return [nearest_pickle_type(SpT) for SpT in spec_type]


def make_h4rg_pixel_map(dead_pix=1, dead_line=1, filename=None):
    """
    Generate a 4096 x 4096 pixel map with a certain percentage of dead pixels
    and lines

    Parameters
    ==========
    dead_pix : int, optional
        The percentage [%] of dead pixels in the frame
    dead_line : int, optional
        The percentage [%] of dead lines in the frame
    filename : str, optional
        Default = None. If None, the array is returned. If not None, a FITS file
        with the pixel map is saved to this path
    """

    chip = np.ones((4096, 4096), dtype=np.float32)
    l = 4096
    w = 4096 // 32

    for i in range(0, 4096, l):
        for j in range(0, 4096, w):
            seg = np.ones((l, w))
            # Add the dead pixels
            x = np.random.randint(l, size=int(w*l*dead_pix/100))
            y = np.random.randint(w, size=int(w*l*dead_pix/100))
            seg[x, y] = 0
            # Add the dead lines
            y = np.random.randint(w, size=int(w*dead_line/100))
            seg[:, y] = 0
            chip[i:i+l, j:j+w] = seg

    if filename is not None:
        hdu = fits.PrimaryHDU(chip)
        hdu.header["SIMCADO"] = ("PIXEL_MAP",
                                 "A map of dead pixels on the detector")
        hdu.writeto(filename, clobber=True)
    else:
        return chip
