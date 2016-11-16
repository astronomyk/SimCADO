"""
simulation.py
"""

# from astropy.io import ascii as ioascii ## unused (OC)
import numpy as np
import simcado as sim

def run(src, mode="wide", cmds=None, opt_train=None, fpa=None,
        detector_layout="small", filename=None, return_internals=False,
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
        ["small", "wide", "zoom", "centre", "full"] Default is "small".
        "small"   - 1x 1k-detector centred in the FoV
        "centre"  - 1x 4k-detector centred in the FoV
        "wide"    - 9x 4k-detector as per MICADO wide field mode (4mas)
        "zoom"    - 9x 4k-detector as per MICADO zoom mode (1.5mas)
        "full"    - "wide" or "zoom" depending on "mode" keyword.

    filename : str, optional
        The filepath for where the FITS images should be saved.
        Default is None. If None, the output images are returned to the user as
        FITS format astropy.io.HDUList objects.

    return_internals : bool
        [False, True] Default is False. If True, the `UserCommands`,
        `OpticalTrain` and `Detector` objects used in the simulation are
        returned in a tuple: `return hdu, (cmds, opt_train, fpa)`

    """

    if cmds is None:
        cmds = sim.UserCommands()
        cmds["INST_FILTER_TC"] = "Ks"
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

    # update any remaining keywords
    cmds.update(kwargs)

    if opt_train is None:
        opt_train = sim.OpticalTrain(cmds)
    if fpa is None:
        fpa = sim.Detector(cmds, small_fov=False)
        
    print(fpa.layout)
    print("Creating", len(cmds.lam_bin_centers), "layer(s) per chip")
    print(len(fpa.chips), "chip(s) will be simulated")
    
    src.apply_optical_train(opt_train, fpa)

    if filename is not None:
        if cmds["OBS_SAVE_ALL_FRAMES"] == "yes":
            for n in cmds["OBS_NDIT"]:
                fname = filename.replace(".",str(n)+".")
                fpa.read_out(filename=fname, to_disk=True, OBS_NDIT=1)
        else:
            fpa.read_out(filename=filename, to_disk=True)
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
    ratio or a list of magnitudes in `mags` in a certain broadband
    `filter_name`.
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
        Number of readouts during the period `exptime`. Default is 1
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