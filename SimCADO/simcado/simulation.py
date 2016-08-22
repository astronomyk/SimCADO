"""
simulation.py
"""

# from astropy.io import ascii as ioascii ## unused (OC)
import numpy as np
import simcado as sim

def run(src, cmds=None, opt_train=None, fpa=None,
        micado_fpa=False, filename=None, return_internals=False,
        **kwargs):
    """
    Run a MICADO simulation with default parameters
    """


    if cmds is None:
        cmds = sim.UserCommands()
        cmds["INST_FILTER_TC"] = "K"
    cmds.update(kwargs)

    if micado_fpa:
        cmds["FPA_CHIP_LAYOUT"] = "default"

    if opt_train is None:
        opt_train = sim.OpticalTrain(cmds)
    if fpa is None:
        fpa = sim.Detector(cmds)

    print(fpa.layout)
    src.apply_optical_train(opt_train, fpa)

    if filename is not None:
        fpa.read_out(filename=filename, to_disk=True)
    else:
        hdu = fpa.read_out()
        if return_internals:
            return hdu, (cmds, opt_train, fpa)
        else:
            return hdu


def snr(mags, filter_name="K", exptime=18000, ndit=1, cmds=None):
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
        Default is "K". Acceptable broadband filters are UBVRIzYJHKKs
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
    cmd = sim.UserCommands()
    cmd["OBS_EXPTIME"] = exptime / ndit
    cmd["OBS_NDIT"] = ndit
    cmd["INST_FILTER_TC"] = filter_name

    opt = sim.OpticalTrain(cmd)

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
