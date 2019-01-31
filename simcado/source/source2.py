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

import numpy as np

from astropy.table import Table
from astropy import units as u
from astropy.io import ascii as ioascii
from astropy.io import fits
from astropy.wcs import WCS

from synphot import SourceSpectrum, SpectralElement
from synphot.models import Empirical1D

from .. import utils


class Source:
    def __init__(self, spectra=None, filename=None, table=None, image=None,
                 lam=None, x=None, y=None, ref=None, weight=None, **kwargs):
        self.meta = {}
        self.meta.update(kwargs)
        
        self.positions = []
        self.spectra = []

        if filename is not None and spectra is not None:
            self._from_file(filename, spectra)

        elif table is not None and spectra is not None:
            self._from_table(table, spectra)

        elif image is not None and spectra is not None:
            self._from_image(image, spectra)

        elif x is not None and y is not None and \
                ref is not None and weight is not None and \
                lam is not None and spectra is not None:
            self._from_arrays(x, y, ref, weight, lam, spectra)

    def _from_file(self, filename, spectra):
        if isinstance(spectra, SourceSpectrum):
            spectra = [spectra]

        if utils.is_fits(filename):
            data = fits.getdata(filename)
            hdr = fits.getheader(filename)
            if isinstance(data, np.ndarray):
                hdu = fits.ImageHDU(data=data, header=hdr)
                hdu.header["ref"] += len(self.spectra)
            else:
                hdr1 = fits.getheader(filename, 1)
                hdr.update(hdr1)
                hdu = fits.BinTableHDU(data=data, header=hdr)
                hdu.data["ref"] += len(self.spectra)










    def _from_table(self, tbl, spectra):
        if isinstance(spectra, SourceSpectrum):
            spectra = [spectra]

        len_spec = len(self.spectra)
        print(len_spec, tbl)
        tbl["ref"] += len_spec
        self.positions += [tbl]
        self.spectra += spectra

    def _from_image(self, image, spectra):
        pass
    
    def _from_arrays(self, x, y, ref, weight, lam, spectra):
        pass

    def image_in_range(self, wave_0, wave_1):
        pass

    def photons_in_range(self, wave_0, wave_1):
        pass

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

    def shift(self, dx, dy):
        pass

    def rotate(self, angle):
        pass

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
