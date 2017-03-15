"""
simulation.py
"""

# from astropy.io import ascii as ioascii ## unused (OC)
import warnings, logging

import numpy as np

import simcado as sim

__all__ = ["run", "snr", "check_chip_positions"]

def run(src, mode="wide", cmds=None, opt_train=None, fpa=None,
        detector_layout="small", filename=None, return_internals=False,
        filter_name=None, exptime=None,
        **kwargs):
    """
    Run a MICADO simulation with default parameters

    Parameters
    ----------
    src : simcado.Source
        The object of interest

    mode : str, optional
        ["wide", "zoom"] Default is "wide", for a 4mas FoV. "Zoom" -> 1.5mas

    cmds : simcado.UserCommands, optional
        A custom set of commands for the simulation. Default is None

    opt_train : simcado.OpticalTrain, optional
        A custom optical train for the simulation. Default is None

    fpa : simcado.Detector, optional
        A custom detector layout for the simulation. Default is None

    detector_layout : str, optional
        ["small", "centre", "full"] Default is "small".
        "small"   - 1x 1k-detector centred in the FoV
        "centre"  - 1x 4k-detector centred in the FoV
        "full"    - "wide" or "zoom" depending on "mode" keyword.

    filename : str, optional
        The filepath for where the FITS images should be saved.
        Default is None. If None, the output images are returned to the user as
        FITS format astropy.io.HDUList objects.

    return_internals : bool
        [False, True] Default is False. If True, the ``UserCommands``,
        ``OpticalTrain`` and ``Detector`` objects used in the simulation are
        returned in a tuple: ``return hdu, (cmds, opt_train, fpa)``

    filter_name : str, TransmissionCurve
        Analogous to passing INST_FILTER_TC as a keyword argument

    exptime : int, float
        [s] Analogous to passing OBS_EXPTIME as a keyword argument

    """

    if cmds is None:
        cmds = sim.UserCommands()
        cmds["INST_FILTER_TC"] = "Ks"

        if detector_layout.lower() in ("small", "centre", "center"):
            cmds["FPA_CHIP_LAYOUT"] = detector_layout

        if mode == "wide":
            cmds["SIM_DETECTOR_PIX_SCALE"] = 0.004
            cmds["INST_NUM_MIRRORS"] = 11
            if detector_layout.lower() == "full":
                cmds["FPA_CHIP_LAYOUT"] = "wide"
        elif mode == "zoom":
            cmds["SIM_DETECTOR_PIX_SCALE"] = 0.0015
            cmds["INST_NUM_MIRRORS"] = 13
            if detector_layout.lower() == "full":
                cmds["FPA_CHIP_LAYOUT"] = "zoom"
        else:
            raise ValueError("'mode' must be either 'wide' or ' zoom', not " + mode)

    if filter_name is not None:
        cmds["INST_FILTER_TC"] = filter_name

    if exptime is not None:
        cmds["OBS_EXPTIME"] = exptime

    # update any remaining keywords
    cmds.update(kwargs)

    if opt_train is None:
        opt_train = sim.OpticalTrain(cmds)
    if fpa is None:
        fpa = sim.Detector(cmds, small_fov=False)

    print("Detector layout")
    print(fpa.layout)
    print("Creating", len(cmds.lam_bin_centers), "layer(s) per chip")
    print(len(fpa.chips), "chip(s) will be simulated")

    src.apply_optical_train(opt_train, fpa)

    if filename is not None:
        if cmds["OBS_SAVE_ALL_FRAMES"] == "yes":
            for n in cmds["OBS_NDIT"]:
                fname = filename.replace(".",str(n)+".")
                hdu = fpa.read_out(filename=fname, to_disk=True, OBS_NDIT=1)
        else:
            hdu = fpa.read_out(filename=filename, to_disk=True)
    else:
        hdu = fpa.read_out()

    if return_internals:
        return hdu, (cmds, opt_train, fpa)
    else:
        return hdu


def snr(mags, filter_name="Ks", total_exptime=18000, ndit=1, cmds=None):
    """
    Return the signal-to-noise for a list of magnitudes in a specific filter

    Uses the standard setup for MICADO and calculates the signal-to-noise
    ratio or a list of magnitudes in ``mags`` in a certain broadband
    ``filter_name``.
    A custom UserCommands object can also be used. Note that this runs a basic
    SimCADO simulation len(mags) times, so execution time can be many minutes.

    Parameters
    ----------
    mags : array-like
        [vega mags] The magnitude(s) of the source(s)

    filter_name : str, optional
        Default is "Ks". Acceptable broadband filters are UBVRIzYJHKKs

    exptime : float
        [s] Total exposure time length. Default is 18000s (5 hours)

    ndit : int, optional
        Number of readouts during the period ``exptime``. Default is 1

    cmds : simcado.UserCommands, optional
        A custom set of commands for the simulations. If not specified, SimCADO
        uses the default MICADO parameters

    Returns
    -------
    sn : np.ndarray
        An array of signal-to-noise ratios for the magnitudes given

    """
    ## TODO: What about argument cmds? (OC)
    if cmds is None:
        cmd = sim.UserCommands()
    else:
        cmd = cmds
    cmd["OBS_EXPTIME"] = total_exptime / ndit
    cmd["OBS_NDIT"] = ndit
    cmd["INST_FILTER_TC"] = filter_name

    opt = sim.OpticalTrain(cmd)

    if type(mags) not in (list, tuple, np.adarray):
        mags = [mags]

    sn = []
    for mag in mags:
        src = sim.source.star(mag)

        fpa = sim.Detector(cmd)
        src.apply_optical_train(opt, fpa, chips=0)
        hdu = fpa.read_out()

        im = hdu[0].data
        cx, cy = np.array(im.shape) // 2
        n = 5
        sig = np.sum(im[cx-n:cx+n+1, cy-n:cy+n+1])
        av = np.average(im[:200, :50])
        # std = np.std(im[:200, :50])    ## unused (OC)

        n_pix = (2*n+1)**2
        only_sig = sig - av*n_pix
        only_noise = av# * np.sqrt(n_pix)  ## TODO: incorrect (OC)

        sn += [only_sig/only_noise]

    return np.array(sn)


def check_chip_positions(filename="src.fits", x_cen=17.084, y_cen=17.084,
                         n = 0.3, mode="wide"):
    """
    Creates a series of grids of stars and generates the output images

    THe number of stars in each grid corresponds to the id number of the chip
    """


    x = [-x_cen]*1 + [0]*2 + [x_cen]*3 + \
        [-x_cen]*4 + [0]*5 + [x_cen]*6 + \
        [-x_cen]*7 + [0]*8 + [x_cen]*9

    y = [-y_cen + i*n for i in range(1)] + \
        [-y_cen + i*n for i in range(2)] + \
        [-y_cen + i*n for i in range(3)] + \
        [0 + i*n for i in range(4)] + \
        [0 + i*n for i in range(5)] + \
        [0 + i*n for i in range(6)] + \
        [y_cen + i*n for i in range(7)] + \
        [y_cen + i*n for i in range(8)] + \
        [y_cen + i*n for i in range(9)]

    lam, spec = sim.source.SED("A0V", "Ks", 15)
    src = sim.source.Source(lam=lam, spectra=spec, x=x, y=y, ref=[0]*len(x))

    sim.run(src, detector_layout="full", filename=filename, mode=mode)

    
    
def _make_snr_grid_fpas(filter_names=["J", "H", "Ks"], mmin=22, mmax=32, 
                        cmds=None, **kwargs):
    """
    Makes a series of :class:`.Detector` objects containing a grid of stars
    
    
    Parameters
    ----------
    filter_names : list
        Which filters to use for the images. See ``simcado.optices.get_filter_set()``
    
    mmin, mmax : float
        [mag] Minimum and maximum magnitudes to use for the grid of stars
        
    cmds : simcado.UserCommands
        A custom set of commands for building the optical train
        
    Optional Parameters
    -------------------
    Any Keyword-Value pairs accepted by a 
    :class:`~simcado.commands.UserCommands` object
    
    Returns
    -------
    fpas : list
        A list of :class:`Detector` objects with the grid of stars for each filter
        len(fpas) == len(filter_names)
    grid : simcado.Source
        A :class:`Source` object containing the grids of stars
    
    See Also
    --------
    :class:`~simcado.commands.UserCommands`
    
    """
    
    if isinstance(filter_names, str):
        filter_names = [filter_names]
    
    if not isinstance(cmds, list):
        cmds = [cmds] * len(filter_names)
        
    fpas = []
    grids = []
    for filt, cmd in zip(filter_names, cmds):
        if cmd is None:
            cmd = sim.UserCommands()
        #cmd["FPA_USE_NOISE"] = "no"
        cmd["OBS_NDIT"] = 1
        cmd["FPA_LINEARITY_CURVE"] = "none"
        cmd["FPA_CHIP_LAYOUT"] = "small"
        cmd.update(kwargs)
        
        grid = sim.source.star_grid(100, mmin, mmax, filter_name=filt, separation=0.4)
        grids += [grid]
        
        hdus, (cmd, opt, fpa) = sim.run(grid,  filter_name=filt, cmds=cmd, return_internals=True)
        fpas += [fpa]
        
    return fpas, grid


def _get_limiting_mags(fpas, grid, exptimes, filter_names=["J", "H", "Ks"], 
                       mmin=22, mmax=32, AB_corrs=None, limiting_sigma=5):
    """
    Return the limiting magnitude(s) for filter(s) and exposure time(s)
    
    
    Parameters
    ----------
    fpas : list
        The output from A list of :class:`Detector` objects with the grid of stars 
        for each filter
    
    grid : simcado.Source
        The :class:`Source` object containing the grid of stars - used for the pixel
        positions of the stars
    
    exptimes : array
        [s] An array of exposure times in seconds
        
    filter_names : list
        A list of filters. See :func:`simcado.optics.get_filter_set`
        
    mmin, mmax : float
        [mag] the minimum and maximum magnitudes in the grid of stars
        
    AB_corrs : list
        [mag] A list of magnitude corrections to convert from Vega to AB magnitudes
        
    limiting_sigma : float
        [\sigma] The number of sigmas to use to define the limiting magnitude. 
        Default is 5*sigma
        
        
    Returns
    -------
    mags_all : list
        [mag] A list of limiting magnitudes for each exposure time for each filter
        Dimensions are [n, m] where n is the number of filters and m is the number
        of exposure times passed
    
    
    """
    
    from scipy import stats
    
    
    if AB_corrs is None:
        AB_corrs = np.zeros(len(fpas))
    if np.isscalar(AB_corrs):
        AB_corrs = [AB_corrs]*len(fpas)
        
        
    if isinstance(filter_names, str):
        filter_names = [filter_names]*len(fpas)
        
    if np.isscalar(exptimes):
        exptimes = [exptimes]*len(fpas)
    
    mags_all = []
    for fpa, filt, AB_corr in zip(fpas, filter_names, AB_corrs):

        lim_mags = []
        for exptime in exptimes:

            hdus = fpa.read_out(OBS_EXPTIME=exptime)
            im = hdus[0].data

            im_width = hdus[0].data.shape[0]
            #x = (grid._x_pix+im_width//2).astype(int)
            x = grid._x_pix.astype(int)
            y = grid._y_pix.astype(int)

            sigs, nss, snrs, bgs = [], [], [], []
            for n in range(len(x)):
                dw = 5
                w = max(dw+5, int((1.-n/len(x))*20))

                ps = im[y[n]-w:y[n]+w, x[n]-w:x[n]+w]

                sig = np.copy(ps[dw:-dw, dw:-dw])
                bg  = np.copy(ps)
                bg[dw:-dw, dw:-dw] = 0

                bgs  += [np.average(bg[bg!=0])]
                nss  += [np.std(bg[bg!=0]) * np.sqrt(np.sum(bg!=0))]
                sigs += [np.sum(sig - bgs[-1])]

            nss  = np.array(nss)
            sigs = np.array(sigs)
            snr = sigs/nss

            mags = np.linspace(mmin, mmax, len(x)) + AB_corr

            mask = snr > 5
            try:
                q = stats.linregress(mags[mask], np.log10(snr[mask]))
                lim_mag = (np.log10(limiting_sigma) - q.intercept) / q.slope
            except:
                lim_mag = 0

            lim_mags += [lim_mag]
            print(exptime, filt, lim_mag)
        mags_all += [lim_mags]
    
    return mags_all


def plot_exptime_vs_limiting_mag(exptimes, limiting_mags, filter_names=["J", "H", "Ks"], 
                                 colors="bgrcymk", mmin=22, mmax=29,
                                 legend_loc=3, marker="+"):
    """
    Plots exposure time versus limiting magnitudes
       
    
    Parameters
    ----------
    exptimes : list, array
        [s] Exposure times corresponding to the signal-to-noise values
    
    limiting_mags : array, list of array
        [mag] Limiting magnitudes for one, or more, filters for the given exposure times
        Dimensions are (1, n) for a single filter, or (m, n) for m filters
        
    filter_names : list
        A list of m filters. See :func:`simcado.optics.get_filter_set`
    
    colors : list
        The colours to use for dots in the plot
        
    mmin, mmax : float
        The minimum and maximum magnitudes for the y axis
        
    marker : str
        The matplotlib scatter marker key
        
    legend_loc : int
        Location of the legend. If ``None`` is passed, no legend is plotted
    
    """

    import matplotlib.pyplot as plt
    
    if len(np.shape(limiting_mags)) == 1:
        limiting_mags = [limiting_mags]
    if filter_names is None:
        filter_names = ["Filter "+str(i) for i in range(np.shape(limiting_mags)[0])]
    
    elif isinstance(filter_names, str):
        filter_names = [filter_names]*np.shape(limiting_mags)[0]

    exptimes = np.array(exptimes)
    
    
    fig = plt.gcf()
    
    #ax = fig.add_axes([a_left, a_bottom, ax_width, ax_height])
    ax1 = fig.add_axes([0, 0, 1, 1])

    for mag, clr, filt in zip(limiting_mags, colors, filter_names):
        plt.plot(exptimes/3600, mag, clr+marker, label=filt)

    if legend_loc is not None:
        plt.legend(loc=legend_loc, scatterpoints=1)
        
    plt.xlabel("Exposure time [hours]")
    plt.ylabel("Limiting Magnitudes")
    plt.xlim(np.min(exptimes/3600)-0.1, np.max(exptimes/3600)+0.1)
    plt.ylim(22,31)

    plt.grid("on")

    ax2 = fig.add_axes([0.5, 0.15, 0.45, 0.35])

    for mag, clr in zip(limiting_mags, colors):
        plt.plot(exptimes, mag, clr+marker)

    plt.plot((60*1,60*1),    (mmin, mmax), "k:")
    plt.text(60*1-5, mmin+0.5, "1 min", horizontalalignment="right")
    plt.plot((60*4, 60*4),   (mmin, mmax), "k:")
    plt.text(60*4-5, mmin+0.5, "4 min", horizontalalignment="right")
    plt.plot((60*15, 60*15), (mmin, mmax), "k:")
    plt.text(60*15-5, mmin+0.5, "15 min", horizontalalignment="right")

    plt.xlim(10, 1800); plt.ylim(mmin, mmax)
    plt.semilogx()
    plt.xlabel("Exposure time [sec]")
    
    

def limiting_mags(exptimes=[1,60,3600,18000], filter_names=["J", "H", "Ks"], 
                  AB_corrs=None, limiting_sigma=5,
                  return_mags=True, make_graph=False, 
                  mmin=22, mmax=31, 
                  cmds=None, **kwargs):
    """
    Return or plot a graph of the limiting magnitudes for MICADO
    
    
    Parameters
    ----------
    exptimes : array
        [s] Exposure times for which limiting magnitudes should be found
    
    filter_names : list
        A list of filters. See :func:`simcado.optics.get_filter_set`
    
    AB_corrs : list
        [mag] A list of magnitude corrections to convert from Vega to AB magnitudes
        
    limiting_sigma : float
        [\sigma] The number of sigmas to use to define the limiting magnitude. 
        Default is 5*sigma
    
    return_mags : bool
        If True (defualt), the limiting magnitude are returned
    
    make_graph : bool
        If True (defualt), a graph of the limiting magnitudes vs exposure time is plotted
        Calls :func:`plot_exptime_vs_limiting_mag`
        
    cmds : simcado.UserCommands
        A custom set of commands for building the optical train    
        
    
    Optional Parameters
    -------------------
    Any Keyword-Value pairs accepted by a :class:`~simcado.UserCommands` object
    
    
    Returns
    -------
    mags_all : list
        [mag] If ``return_mags=True``, returns a list of limiting magnitudes for 
        each exposure time for each filter
        Dimensions are [n, m] where n is the number of filters and m is the number
        of exposure times passed
    
    
    Notes
    -----
    Vega to AB = {"J" : 0.91 , "H" : 1.39 , "Ks" : 1.85}
    
    
    Examples
    --------
    :
        >>> # Set 30 logarithmic time bins between 1 sec and 5 hours
        >>> exptimes = np.logspace(0, np.log10(18000), num=30, endpoint=True)
        >>> limiting_mags(exptimes=exptimes, filter_names=["J", "PaBeta"], 
        ...               make_graph=False)
    
    """


    fpas, grid    = _make_snr_grid_fpas(filter_names, cmds=cmds, 
                                        mmin=mmin, mmax=mmax, **kwargs)    
    limiting_mags = _get_limiting_mags(fpas, grid, exptimes, filter_names, 
                                       mmin=mmin, mmax=mmax, AB_corrs=AB_corrs,
                                       limiting_sigma=limiting_sigma)
    
    if make_graph:
        plot_exptime_vs_limiting_mag(exptimes, limiting_mags, filter_names, 
                                     mmin=mmin, mmax=mmax)
    
    if return_mags:
        return limiting_mags    