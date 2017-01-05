"""
A module that holds the functions needed to plot the increase in flux due to
non-zero blocking outside of a filter's wavelength range

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import simcado as sim



def plot_filter_wings_flux(lam_min=1.490 , lam_max=1.780, spec_types=["A0V", "G2V", "M5V"],
                           plot_atmosphere=True, **kwargs):
    """
    Plot the effect of non-zero blocking in the filter wings
    
    Parameters
    ----------
    lam_min, lam_max : float
        [um] the cut-off wavelengths for the filter
    spec_type : list of strings, optional
        The spectral types to be used for the comparison
    plot_atmosphere : bool
        Default is True. Plot only the atmospheric emission, i.e. no star
    
    Optional parameters
    -------------------
    Any keyword-value pair from a simcado config file
    
    Notes
    -----
    The current filter fluxes are calculated with all sources of noise turned off for 
    the stars. No shot noise is taken into accoung either.
    For the atmosphere, all background emission is turned on, but still no shot noise is used
    """
    
    kwargs = {"ATMO_USE_ATMO_BG "     : "no", 
              "SCOPE_USE_MIRROR_BG"   : "no",
              "INST_USE_AO_MIRROR_BG" : "no", 
              "FPA_USE_NOISE"         : "no"}
    
    filts = make_tophat_tcs(lam_min, lam_max)

    
    src_stars = [sim.source.star(20, "H", spec_type=spt) for spt in spec_types]
    stars = get_flux(src_stars, filts, **kwargs)

    plt.figure(figsize=(10,10))
    plt.suptitle("Narrow band H-cont filter (1.5685um to 1.5915um)", fontsize=20)

    for i, ttl in zip(range(len(spec_types)), spec_types):

        star_tbl = make_wing_tbl(stars[i])

        plt.subplot(2, len(spec_types)%2+1, i+1)
        plot_flux_vs_wing(star_tbl, loc=4)
        plt.title(ttl)
        plt.ylim(1E-5,1E-1)
        
    if plot_atmosphere:
        
        src_sky = [sim.source.empty_sky()]
        sky = get_flux(src_sky, filts)
        
        sky_tbl = make_wing_tbl(sky)
        
        plt.subplot(2, len(spec_types)%2+1, len(spec_types)+1)
        plot_flux_vs_wing(sky_tbl, loc=4)
        plt.title("Atmospheric BG")
        plt.ylim(1E-5,1E-1)


def get_flux(srcs, filts, **kwargs):
    """
    Returns the total flux from a series of sources through a series of filters 
    
    Parameters
    ----------
    srcs : list of simcado.Source objects
    filts : list of simcado.TransmissionCurve objects
    
    Optional Parameters
    -------------------
    kwargs : any keyword-value pair for a simcado.UserCommands object
    
    Returns
    -------
    fluxes : list of floats
        A list of fluxes per filter for each source
    
    """
    
    cmd = sim.UserCommands()
    cmd.update(kwargs)
    
    fluxes = []
    for src in srcs:
        flux_i = []
        for filt in filts:
            
            cmd["INST_FILTER_TC"] = filt
            opt = sim.OpticalTrain(cmd)
            fpa = sim.Detector(cmd)

            src.apply_optical_train(opt, fpa)
            flux_i += [np.sum(fpa.chips[0].array)]

        fluxes += [flux_i]
   
    if len(fluxes) == 1:
        fluxes = fluxes[0]
    return fluxes


def make_tophat_tcs(lam_min, lam_max):
    """
    Create three simcado.TransmissionCurve objects for a filter and the red/blue wings
    
    Parameters
    ----------
    lam_min, lam_max : float
        The blue and red cutoffs for the filter
        
    Returns
    -------
    tc_blue, tc_filt, tc_red : simcado.TransmissionCurves
        FIlter curves for the filter, plus "filter" curves for the red and blue wings
    
    """

    lam  = np.array([0.3, lam_min-0.001, lam_min, 
                          lam_max, lam_max+0.001, 3.])
    blue = np.array([1,1,0,0,0,0])
    filt = np.array([0,0,1,1,0,0])
    red  = np.array([0,0,0,0,1,1])

    tc_blue = sim.spectral.TransmissionCurve(lam=lam, val=blue, lam_res=0.001)
    tc_filt = sim.spectral.TransmissionCurve(lam=lam, val=filt, lam_res=0.001)
    tc_red  = sim.spectral.TransmissionCurve(lam=lam, val=red , lam_res=0.001)
    
    return tc_blue, tc_filt, tc_red


def make_wing_tbl(flux, trans = [1E-2, 1E-3, 1E-4, 1E-5]):
    """
    Make an astropy.Table with the relative flux increase per wing for a series of blocking coeffients
    
    Flux is a list of 3 values: the first and third are the fluxed through the blue and red wings, the second
    middle value is the flux through the filter wavelength range. All fluxes should be for a 100% transmission.
    The function iterates over the values in `trans`  (transmission coefficients) to make a table with the
    relative fluxes between filter and the (blue/both/red) wing for a certain blocking coefficient.
    
    Parameters
    ----------
    flux : list of 3 floats
        list with the fluxes for the filter and the two wings
    trans : list of floats, optional
        the blocking coefficients for the wings
        
    Returns
    -------
    wing_tbl : astropy.Table
    
    """
    
    from astropy.table import Table
    
    blue = [f*flux[0]/flux[1] for f in trans]
    red  = [f*flux[2]/flux[1] for f in trans]
    both = [f*(flux[0]+flux[2])/flux[1] for f in trans]

    wing_tbl = Table(data=(trans, blue, both, red), names=["Wing_Factor", "Blue_Wing", "Both_Wings","Red_Wing"])
    return wing_tbl


def plot_flux_vs_wing(tbl, loc=2):
    """
    Plots the relative flux increase due to different wing blocking coefficients.
    
    Plots to the current axis. Subplot should be set beforehand and show() called afterwards
    
    Parameters
    ----------
    tbl : astropy.Table
        The output from make_wing_tbl
    loc : int, optional
        location for the legend
        
    """
    for col, clr in zip(tbl.colnames[1:], "bkr"):
        plt.plot(tbl["Wing_Factor"], tbl[col], clr, label=col)
        plt.loglog()
        plt.grid("on")
        plt.legend(loc=loc)
        plt.xlabel("Wing transmission")
        plt.ylabel("Fractional flux increase")