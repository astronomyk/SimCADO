import glob
import os

import numpy as np
from astropy.io import fits
import poppy

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def display_profiles(HDUlist_or_filename=None, ext=0, overplot=False,
                     title=None, vmin=1e-8, vmax=1e-2, **kwargs):
    """ Produce two plots of PSF radial profile and encircled energy

    !!! Taken directly from Poppy !!!

    See also the display_ee function.

    Parameters
    ----------
    HDUlist_or_filename1,2 : fits.HDUlist or string
        FITS files containing image to difference
    ext : bool
        FITS extension to use. Default is 0
    overplot : bool
        whether to overplot or clear and produce an new plot. Default false
    title : string, optional
        Title for plot

    """
    if isinstance(HDUlist_or_filename, str):
        hdu_list = fits.open(HDUlist_or_filename, ext=ext)
    elif isinstance(HDUlist_or_filename, fits.HDUList):
        hdu_list = HDUlist_or_filename
    else:
        raise ValueError("input must be a filename or HDUlist")

    radius, profile, ee = poppy.radial_profile(hdu_list, ee=True, ext=ext,
                                               **kwargs)
    plt.semilogy(radius, profile)

    plt.xlabel("Radius [arcsec]")
    plt.ylabel("PSF radial profile")
    plt.xlim(0, 1)
    plt.ylim(vmin, vmax)

    plt.twinx()

    plt.plot(radius, ee, color='r')

    # fwhm = 2 * radius[np.where(profile < profile[0] * 0.5)[0][0]]
    # plt.text(fwhm, 0.9, 'FWHM = %.3f"' % fwhm)

    plt.xlabel("Radius [arcsec]")
    plt.ylabel("Encircled Energy")

    for level in [0.5, 0.8, 0.95]:
        if (ee > level).any():
            ee_lev = radius[np.where(ee > level)[0][0]]
            yoffset = 0 if level < 0.9 else -0.05
            plt.text(ee_lev + 0.1, level + yoffset,
                     'EE=%2d%% at r=%.3f"' % (level * 100, ee_lev))


def plot_psf_summary(filename, ext, layer, **kwargs):

    with fits.open(filename) as hdulist:
        hdu = hdulist[ext]
        if "CATTYPE" in hdu.header:
            return

        if len(hdu.data.shape) == 3:
            data = hdu.data[layer, :, :]
        else:
            data = hdu.data

        data /= np.sum(data)

        prihdu = fits.PrimaryHDU(header=hdu.header, data=data)
        if "deg" in prihdu.header["CUNIT1"]:
            prihdu.header["CDELT1"] *= 3600
            prihdu.header["CDELT2"] *= 3600
            prihdu.header["CUNIT1"] = "arcsec"
            prihdu.header["CUNIT2"] = "arcsec"
        prihdu.header["PIXELSCL"] = prihdu.header["CDELT1"]

        hdulist = fits.HDUList([prihdu])

    wave = hdulist[0].header["WAVE0"]
    title = "{}[{}]({}) - {} um" \
            "".format(os.path.basename(psf_file), ext, layer, wave)
    print(title)

    params = {"HDUlist_or_filename": hdulist,
              "pixelscale": "CDELT1",
              "interpolation": "none",
              "normalize": "total",
              "title": "",
              "vmax": 1e-2,
              "vmin": 1e-8,
              "cmap": "viridis"}
    params.update(kwargs)

    plt.figure(figsize=(10, 8))
    plt.suptitle(title)

    plt.axes([0.05, 0.55, 0.4, 0.4])
    poppy.display_psf(**params, imagecrop=2)
    plt.xlabel("[arcsec]")
    plt.ylabel("[arcsec]")

    plt.axes([0.55, 0.55, 0.4, 0.4])
    poppy.display_psf(**params, imagecrop=0.2)
    plt.xlabel("[arcsec]")
    plt.ylabel("[arcsec]")

    plt.axes([0.05, 0.05, 0.9, 0.4])
    display_profiles(hdulist, title=title, vmin=params["vmin"],
                     vmax=params["vmax"])

    return title


psf_list = glob.glob("C:/Work/irdb/_PSFs/*.fits")

for psf_file in psf_list:
    info = fits.info(psf_file, output=False)
    for ii, hdu_info in enumerate(info):
        if hdu_info[3] == "BinTableHDU":
            continue
        if len(hdu_info[5]) >= 2:
            ext = ii

            if "MAORY_SCAO_FVPSF_1.5mas" in psf_file:
                layer = 4
            elif "MAORY_SCAO_FVPSF_4mas" in psf_file:
                layer = 24
            else:
                layer = 0

            title = plot_psf_summary(psf_file, ext, layer)
            if title is None:
                continue
            plt.savefig("C:/Work/irdb/_PSFs/images/psf_summary/{}.png"
                        "".format(title))
