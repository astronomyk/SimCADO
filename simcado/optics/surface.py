import os
import warnings

import numpy as np

from astropy import units as u
from astropy.io import ascii as ioascii
from astropy.table import Table

from synphot import SpectralElement, SourceSpectrum
from synphot.models import Empirical1D, BlackBody1D
from synphot.units import PHOTLAM

from .. import utils


class SpectralSurface:
    def __init__(self, filename=None, **kwargs):
        self.meta = {"filename"   : filename,
                     "area"       : 0,      # m2
                     "etendue"    : 0,      # arcsec2 * m2
                     "temp"       : -270,   # deg C
                     "emission_unit" : "",
                     "wavelength_unit" : "um"}
        self.meta.update(kwargs)

        self.table = Table()
        if filename is not None and os.path.exists(filename):
            self.table = ioascii.read(filename)
            tbl_meta = utils.convert_table_comments_to_dict(self.table)
            if isinstance(tbl_meta, dict):
                self.meta.update(tbl_meta)

    @property
    def wavelength(self):
        return self._get_array("wavelength")

    @property
    def transmission(self):
        return self._get_ter_property("transmission")

    @property
    def emissivity(self):
        return self._get_ter_property("emissivity")

    @property
    def reflection(self):
        return self._get_ter_property("reflection")

    @property
    def emission(self):

        flux = self._get_array("emission")
        if flux is not None:
            wave = self._get_array("wavelength")
            flux = make_emission_from_array(flux, wave)
        elif "temperature" in self.meta:
            emiss = self.emissivity                # SpectralElement
            temp = self.meta["temperature"] + 273  # BlackBody1D --> Kelvin
            flux = make_emission_from_emissivity(temp, emiss)
        else:
            flux = None

        return flux

    def _get_ter_property(self, ter_property):
        compliment_names = ["transmission", "emissivity", "reflection"]
        ii = np.where(ter_property == np.array([compliment_names]))[0][0]
        compliment_names.pop(ii)

        wave = self._get_array("wavelength")
        value_arr = self._get_array(ter_property)
        if value_arr is None:
            value_arr = self._compliment_array(*compliment_names)
        if value_arr is not None and wave is not None:
            response_curve = SpectralElement(Empirical1D, points=wave,
                                             lookup_table=value_arr)
        else:
            response_curve = None
            warnings.warn("Both wavelength and {} must be set"
                          "".format(ter_property))

        return response_curve

    def _compliment_array(self, colname_a, colname_b):
        col_a = self._get_array(colname_a)
        col_b = self._get_array(colname_b)

        if col_a is not None and col_b is not None:
            col_c = 1*col_a.unit - (col_a + col_b)
        elif col_a is not None and col_b is None:
            col_c = 1*col_a.unit - col_a
        elif col_b is not None and col_a is None:
            col_c = 1*col_b.unit - col_b
        elif col_b is None and col_a is None:
            col_c = None

        return col_c

    def _get_array(self, colname):
        if colname in self.meta:
            val = self.meta[colname]
        elif colname in self.table.colnames:
            val = self.table[colname].data
        else:
            warnings.warn("{} not found in either '.meta' or '.table'"
                          "".format(colname))
            return None

        col_units = colname+"_unit"
        if col_units in self.meta:
            units = u.Unit(self.meta[col_units])
        elif isinstance(val, u.Quantity):
            units = val.unit
        else:
            units = u.Unit("")

        if isinstance(val, u.Quantity):
            val_out = val.to(units)
        elif isinstance(val, (list, tuple, np.ndarray)):
            val_out = val * units
        elif val is None:
            val_out = val
        else:
            raise ValueError("{} must be of type: Quantity, array, list, tuple"
                             "".format(colname))

        return val_out


def quantify(item, unit):
    if isinstance(item, u.Quantity):
        quant = item.to(u.Unit(unit))
    else:
        quant = item * u.Unit(unit)
    return quant


def make_emission_from_emissivity(temp, emiss_src_spec):
    if emiss_src_spec is None:
        warnings.warn("Either emission or emissivity must be set")
        flux = None
    else:
        flux = SourceSpectrum(BlackBody1D, temp) / u.sr
        flux = flux * emiss_src_spec

    return flux


def make_emission_from_array(flux, wave, meta):
    if "emission_unit" not in meta and \
            not isinstance(flux, u.Quantity):
        warnings.warn("emission_unit must be set, or emission must"
                      "be an astropy.Quantity")
        flux = None
    elif "emission_unit" in meta:
        flux = quantify(flux, meta["emission_unit"])

    if isinstance(wave, u.Quantity) and isinstance(flux, u.Quantity):
        flux_unit, angle = extract_type_from_unit(flux.unit, "angle")
        flux = flux / angle

        if is_flux_binned(flux.unit):
            flux = normalise_binned_flux(flux, wave)

        flux = SourceSpectrum(Empirical1D, points=wave,
                              lookup_table=flux)
        flux.meta["angle_unit"] = angle
    else:
        warnings.warn("wavelength and emission must be "
                      "astropy.Quantity objects")
        flux = None
    return flux


def extract_type_from_unit(unit, unit_type):
    unit = unit**1
    extracted_units = u.Unit("")
    for base, power in zip(unit._bases, unit._powers):
        if unit_type in base.physical_type:
            extracted_units *= base**power

    new_unit = unit / extracted_units

    return new_unit, extracted_units


def extract_base_from_unit(unit, base_unit):
    unit = unit**1
    extracted_units = u.Unit("")
    for base, power in zip(unit._bases, unit._powers):
        if base == base_unit:
            extracted_units *= base**power

    return unit * extracted_units**-1, extracted_units


def normalise_binned_flux(flux, wave):

    bins = np.zeros(len(wave)) * wave.unit
    bins[:-1] = 0.5 * np.diff(wave)
    bins[1:] += 0.5 * np.diff(wave)
    # bins[0] *= 2.   # edge bins only have half the flux of other bins
    # bins[-1] *= 2.

    bin_unit = extract_base_from_unit(flux.unit, u.bin)[1]
    flux = flux / bins / bin_unit

    return flux


def is_flux_binned(unit):
    unit = unit**1
    flag = False
    if u.bin in unit._bases or "flux density" not in unit.physical_type:
        flag = True

    return flag





