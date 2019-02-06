# old functionality to implement:
# - provide x, y, lam, spectra, weight, ref
# - overridden + : number, Source, SourceSpectrum
# - overridden * : number, SpectralElement
# - write to and read from file
# - shift all positions
# - rotate around the centre
# - photons_in_range returns the photons per spectrum in a wavelength range
# - image_in_range returns an image of the source for a wavelength range
#
# old functionality which will be removed:
# - project_onto_chip
# - apply_optical_train
#
# old structure --> new structure:
# - all data held in 6 arrays
# --> new dicts for positions, spectrum
#       position can be a Table or an ImageHDU
#       spectrum is a SourceSpectrum
#
# Use cases:
# image + spectrum
# images + spectra
# table + spectrum
# table + spectra
#
# table columns = x, y, spec_id, weight
# table meta keywords = x_unit, y_unit
#
# image header keywords = WCS, SPEC_ID, WEIGHT
# [WCS = CRPIXn, CRVALn = (0,0), CTYPEn, CDn_m, NAXISn, CUNITn

import pickle
import warnings

import numpy as np

from astropy.table import Table, Column
from astropy.io import ascii as ioascii
from astropy.io import fits
from astropy import wcs
from astropy import units as u

from synphot import SourceSpectrum, SpectralElement, Observation
from synphot.models import Empirical1D

from ..optics.image_plane import ImagePlane, make_image_plane_header
from .. import utils


class Source:

    def __init__(self, filename=None,
                 lam=None, spectra=None, x=None, y=None, ref=None, weight=None,
                 table=None, image=None, **kwargs):

        self.meta = {}
        self.meta.update(kwargs)
        
        self.positions = []
        self.spectra = []

        self.bandpass = None

        valid = validate_source_input(lam=lam, x=x, y=y, ref=ref, weight=weight,
                                      spectra=spectra, table=table, image=image,
                                      filename=filename)

        spectra = convert_to_list_of_spectra(spectra, lam)

        if filename is not None and spectra is not None:
            self._from_file(filename, spectra)

        elif table is not None and spectra is not None:
            self._from_table(table, spectra)

        elif image is not None and spectra is not None:
            self._from_image(image, spectra)

        elif x is not None and y is not None and \
                ref is not None and spectra is not None:
            self._from_arrays(x, y, ref, weight, spectra)

    def _from_file(self, filename, spectra):
        filename = utils.find_file(filename)

        if utils.is_fits(filename):
            fits_type = utils.get_fits_type(filename)
            data = fits.getdata(filename)
            hdr = fits.getheader(filename)
            if fits_type == "image":
                image = fits.ImageHDU(data=data, header=hdr)
                self._from_image(image, spectra)
            elif fits_type == "bintable":
                hdr1 = fits.getheader(filename, 1)
                hdr.update(hdr1)
                tbl = Table(data, meta=dict(hdr))
                tbl.meta.update(utils.convert_table_comments_to_dict(tbl))
                self._from_table(tbl, spectra)
        else:
            tbl = ioascii.read(filename)
            tbl.meta.update(utils.convert_table_comments_to_dict(tbl))
            self._from_table(tbl, spectra)

    def _from_table(self, tbl, spectra):
        if "weight" not in tbl.colnames:
            tbl.add_column(Column(name="weight", data=np.ones(len(tbl))))
        tbl["ref"] += len(self.spectra)
        self.positions += [tbl]
        self.spectra += spectra

    def _from_image(self, image, spectra):
        if spectra is not None:
            image.header["SPEC_REF"] = len(self.spectra)
            self.positions += [image]
            self.spectra += spectra
        else:
            image.header["SPEC_REF"] = ""

    def _from_arrays(self, x, y, ref, weight, spectra):
        if weight is None:
            weight = np.ones(len(x))

        x = utils.quantify(x, u.arcsec)
        y = utils.quantify(y, u.arcsec)
        tbl = Table(names=["x", "y", "ref", "weight"],
                    data=[x, y, np.array(ref) + len(self.spectra), weight])
        tbl.meta["x_unit"] = "arcsec"
        tbl.meta["y_unit"] = "arcsec"

        self.positions += [tbl]
        self.spectra += spectra

    def image_in_range(self, wave_min, wave_max, pixel_scale=1*u.arcsec,
                       layers=None, area=None, nbins=100):
        if layers is None:
            layers = list(np.arange(len(self.positions)))

        positions = [self.positions[ii] for ii in layers]
        hdr = make_image_plane_header(positions, pixel_scale=pixel_scale)
        im_plane = ImagePlane(hdr)

        for position in positions:
            if isinstance(position, Table):
                fluxes = self.photons_in_range(wave_min, wave_max,
                                               position["ref"], area, nbins)
                tbl = Table(names=["x", "y", "flux"],
                            data=[position["x"], position["y"], fluxes])
                tbl.meta.update(position.meta)
                hdu_or_table = tbl

            elif isinstance(position, fits.ImageHDU):
                if position.header["SPEC_REF"] is not "":
                    flux = self.photons_in_range(wave_min, wave_max,
                                                 position["ref"], area, nbins)
                else:
                    flux = 1

                image = position.data * flux
                hdu = fits.ImageHDU(header=position.header, data=image)
                hdu_or_table = hdu
            else:
                continue

            im_plane.add(hdu_or_table)

        return im_plane

    def photons_in_range(self, wave_min, wave_max, indexes=None, area=None,
                         nbins=100):
        if indexes is None:
            indexes = range(len(self.spectra))

        spectra = [self.spectra[ii] for ii in indexes]
        counts = photons_in_range(spectra, wave_min, wave_max, area=area,
                                  bandpass=self.bandpass, nbins=nbins)
        return counts

    @classmethod
    def load(cls, filename):
        """Load :class:'.Source' object from filename"""
        with open(filename, 'rb') as fp1:
            src = pickle.load(fp1)
        return src

    def dump(self, filename):
        """Save to filename as a pickle"""
        with open(filename, 'wb') as fp1:
            pickle.dump(self, fp1)

    def write_to_fits(self, filename):
        pass

    def read_from_fits(self, filename):
        pass

    def shift(self, dx=0, dy=0, layers=None):
        """
        Shifts the coordinates of the source by (dx, dy) in [arcsec]

        Parameters
        ----------
        dx, dy : float, array
            [arcsec] The offsets for each coordinate in the arrays ``x``, ``y``.
            - If dx, dy are floats, the same offset is applied to all coordinates
            - If dx, dy are arrays, they must be the same length as ``x``, ``y``

        """

        if layers is None:
            layers = np.arange(len(self.positions))

        for ii in layers:
            if isinstance(self.positions[ii], Table):
                x = utils.quantity_from_table("x", u.arcsec)
                x += utils.quantify(dx, u.arcsec)
                self.positions[ii]["x"] = x

                y = utils.quantity_from_table("y", u.arcsec)
                y += utils.quantify(dy, u.arcsec)
                self.positions[ii]["y"] = y
            elif isinstance(self.positions[ii], fits.ImageHDU):
                dx = utils.quantify(dx, "arcsec").to(u.deg)
                dy = utils.quantify(dy, "arcsec").to(u.deg)
                self.positions[ii].header["CRVAL1"] += dx.value
                self.positions[ii].header["CRVAL2"] += dy.value

    def rotate(self, angle, offset=None, layer=None):
        pass

    def add_bandpass(self, bandpass):
        if not isinstance(bandpass, SpectralElement):
            raise ValueError("type(bandpass) must be synphot.SpectralElement")

        self.bandpass = bandpass

    def __add__(self, other):
        pass

    def __mul__(self, other):
        pass

    def __sub__(self, other):
        pass

    def __radd__(self, other):
        pass

    def __rmul__(self, other):
        pass

    def __rsub__(self, other):
        pass


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





