import warnings

import numpy as np
from astropy import wcs, units as u
from astropy.io import fits
from astropy.table import Table
from synphot import SourceSpectrum, Empirical1D, SpectralElement, Observation

from simcado import utils


def validate_source_input(**kwargs):
    if "filename" in kwargs and kwargs["filename"] is not None:
        filename = kwargs["filename"]
        if utils.find_file(filename) is None:
            warnings.warn("filename was not found: {}".format(filename))

    if "image" in kwargs and kwargs["image"] is not None:
        image = kwargs["image"]
        if not isinstance(image, (fits.PrimaryHDU, fits.ImageHDU)):
            raise ValueError("image must be fits.HDU object with a WCS."
                             "type(image) == {}".format(type(image)))

        if len(wcs.find_all_wcs(image.header)) == 0:
            warnings.warn("image does not contain valid WCS. {}"
                          "".format(wcs.WCS(image)))

    if "table" in kwargs and kwargs["table"] is not None:
        tbl = kwargs["table"]
        if not isinstance(tbl, Table):
            raise ValueError("table must be an astropy.Table object:"
                             "{}".format(type(tbl)))

        if not np.all([col in tbl.colnames for col in ["x", "y", "ref"]]):
            raise ValueError("table must contain at least column names: "
                             "'x, y, ref': {}".format(tbl.colnames))

    return True


def convert_to_list_of_spectra(spectra, lam):
    spectra_list = []
    if isinstance(spectra, SourceSpectrum):
        spectra_list += [spectra]

    elif isinstance(spectra, (tuple, list)) and \
            isinstance(spectra[0], SourceSpectrum):
        spectra_list += spectra

    elif isinstance(spectra, np.ndarray) and isinstance(lam, np.ndarray) and \
            len(spectra.shape) == 1 :
        spec = SourceSpectrum(Empirical1D, points=lam, lookup_table=spectra)
        spectra_list += [spec]

    elif ((isinstance(spectra, np.ndarray) and
           len(spectra.shape) == 2) or
          (isinstance(spectra, (list, tuple)) and
           isinstance(spectra[0], np.ndarray))) and \
            isinstance(lam, np.ndarray):

        for sp in spectra:
            spec = SourceSpectrum(Empirical1D, points=lam, lookup_table=sp)
            spectra_list += [spec]

    return spectra_list


def photons_in_range(spectra, wave_min, wave_max,
                     bandpass=None, area=None, nbins=100):
    wave_min = utils.quantify(wave_min, u.um)
    wave_max = utils.quantify(wave_max, u.um)
    wave = np.linspace(wave_min.value, wave_max.value, nbins) * u.um

    per_unit_area = False
    if area is None:
        area = 1 * u.m ** 2
        per_unit_area = True

    if bandpass is not None:
        bandpass = bandpass
    else:
        bandpass = SpectralElement(Empirical1D, points=wave,
                                   lookup_table=[1] * len(wave))
    counts = [Observation(spectrum,
                          bandpass).effstim(flux_unit="count", binned=True,
                                            area=area, waverange=wave).value
              for spectrum in spectra]
    # effstim doesn't account for the edges well. Generally one too many bins
    correction_factor = len(wave) / (len(wave) + 1.)
    counts *= u.ct / u.s * correction_factor

    if per_unit_area:
        counts = counts / area

    return counts


def make_imagehdu_from_table(x, y, flux, pix_scale=1*u.arcsec):

    pix_scale = pix_scale.to(u.deg)
    unit = pix_scale.unit
    x = utils.quantify(x, unit)
    y = utils.quantify(y, unit)

    xpixmin = int(np.floor(np.min(x) / pix_scale))
    ypixmin = int(np.floor(np.min(y) / pix_scale))
    xvalmin = (xpixmin * pix_scale).value
    yvalmin = (ypixmin * pix_scale).value

    the_wcs = wcs.WCS(naxis=2)
    the_wcs.wcs.crpix = [0., 0.]
    the_wcs.wcs.cdelt = [pix_scale.value, pix_scale.value]
    the_wcs.wcs.crval = [xvalmin, yvalmin]
    the_wcs.wcs.cunit = [unit, unit]
    the_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    ypix, xpix = the_wcs.wcs_world2pix(y.to(u.deg), x.to(u.deg), 1)
    yint, xint  = ypix.astype(int), xpix.astype(int)

    image = np.zeros((np.max(xint) + 1, np.max(yint) + 1))
    for ii in range(len(xint)):
        image[xint[ii], yint[ii]] += flux[ii]

    hdu = fits.ImageHDU(data=image)
    hdu.header.extend(the_wcs.to_header())

    return hdu


def scale_imagehdu(imagehdu, area=None, solid_angle=None, waverange=None):

    if "BUNIT" in imagehdu.header:
        unit = u.Unit(imagehdu.header["BUNIT"])
    elif "FLUXUNIT" in imagehdu.header:
        unit = u.Unit(imagehdu.header["BUNIT"])
    else:
        unit = ""

    zero  = 0  * u.Unit(unit)
    scale = 1 * u.Unit(unit)
    if area is not None:
        scale *= area
    if solid_angle is not None:
        scale *= solid_angle
    if waverange is not None:
        scale *= waverange
    if "BSCALE" in imagehdu.header:
        scale *= imagehdu.header["BSCALE"]
        imagehdu.header["BSCALE"] = 1
    if "BZERO" in imagehdu.header:
        zero = imagehdu.header["BZERO"]
        imagehdu.header["BZERO"] = 0

    imagehdu.data = imagehdu * scale + zero
    imagehdu.header["BUNIT"] = str(imagehdu.data.unit)
    imagehdu.header["FLUXUNIT"] = str(imagehdu.data.unit)

    return imagehdu