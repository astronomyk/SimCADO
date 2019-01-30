import os
import warnings

import numpy as np

from astropy import units as u
from astropy.io import ascii as ioascii
from astropy.table import Table

from synphot import SpectralElement, SourceSpectrum
from synphot.models import Empirical1D, BlackBody1D

from .. import utils


class SpectralSurface:
    """
    Initialised by a file containing one or more of the following columns:
    transmission, emissivity, reflection. The column wavelength must be given.
    Alternatively kwargs for the above mentioned quantities can be passed as
    arrays. If they aren't Quantities, then a unit should also be passed with
    the <array_name>_unit syntax (i.e. emission_unit or wavelength_unit)

    """
    def __init__(self, filename=None, **kwargs):
        filename = utils.find_file(filename)
        self.meta = {"filename"   : filename,
                     "temp"       : -270*u.deg_C,  # deg C
                     "emission_unit" : "",
                     "wavelength_unit" : u.um}

        self.table = Table()
        if filename is not None and os.path.exists(filename):
            self.table = ioascii.read(filename)
            tbl_meta = utils.convert_table_comments_to_dict(self.table)
            if isinstance(tbl_meta, dict):
                self.meta.update(tbl_meta)

        self.meta.update(kwargs)

    @property
    def area(self):
        if "area" in self.meta:
            the_area = self.from_meta("area", u.m**2)
        elif "outer" in self.meta:
            outer_diameter = self.from_meta("outer", u.m)
            the_area = np.pi * (0.5 * outer_diameter)**2
            if "inner" in self.meta:
                inner_diameter = self.from_meta("inner", u.m)
                the_area -= np.pi * (0.5 * inner_diameter) ** 2
        else:
            the_area = None
        return the_area

    @property
    def mirror_angle(self):
        if "angle" in self.meta:
            mirr_angle = self.from_meta("angle", u.deg)
        else:
            mirr_angle = 0
        return mirr_angle

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
        """
        Looks for an emission array in self.meta. If it doesn't find this, it
        defaults to creating a blackbody and multiplies this by the emissivity.
        Assumption is that self.meta["temp"] is in deg_C
        Return units are in PHOTLAM arcsec^-2, even though arcsec^-2 is not
        given
        """

        flux = self._get_array("emission")
        if flux is not None:
            wave = self._get_array("wavelength")
            flux = make_emission_from_array(flux, wave, meta=self.meta)
        elif "temp" in self.meta:
            emiss = self.emissivity                     # SpectralElement [0..1]
            temp = quantify(self.meta["temp"], u.deg_C).value + 273.
            flux = make_emission_from_emissivity(temp, emiss)
        else:
            flux = None

        has_solid_angle = extract_type_from_unit(flux.meta["solid_angle"],
                                                 "solid angle")[1] != u.Unit("")
        if flux is not None and has_solid_angle:
            conversion_factor = flux.meta["solid_angle"].to(u.arcsec ** -2)
            flux = flux * conversion_factor
            flux.meta["solid_angle"] = u.arcsec**-2
            flux.meta["history"] += ["Converted to arcsec-2: {}"
                                     "".format(conversion_factor)]

        return flux

    def from_meta(self, key, default_unit=None):
        """
        Converts a specific value in the meta dict to a ``Quantity``

        Parameters
        ----------
        key : str
            Which key to pull from self.meta
        default_unit : str, Unit
            In case self.meta doesn't contain a unit for the desired key

        Returns
        -------
        meta_quantity : Quantity

        """

        if default_unit is None:
            default_unit = ""
        meta_quantity = get_meta_quantity(self.meta, key, u.Unit(default_unit))

        return meta_quantity

    def _get_ter_property(self, ter_property):
        """
        Looks for arrays for transmission, emissivity, or reflection

        Parameters
        ----------
        ter_property : str
            ``transmission``, ``emissivity``, ``reflection``

        Returns
        -------
            response_curve : ``synphot.SpectralElement``

        """

        compliment_names = ["transmission", "emissivity", "reflection"]
        ii = np.where([ter_property == name for name in compliment_names])[0][0]
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
        """
        Returns an complimentary array using: ``a + b + c = 1``

        E.g. ``Emissivity = 1 - (Transmission + Reflection)``

        Parameters
        ----------
        colname_a : str
            Name of the first TER property to look for
        colname_b
            Name of the second TER property to look for

        Returns
        -------
        col_c : ``synphot.SpectralElement``
            Complimentary spectrum to those given

        """

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
        """
        Looks for an array in either the self.meta or self.table attributes

        Order of search goes: 1. self.meta, 2. self.table

        Parameters
        ----------
        colname : str
            Array column (or key) name

        Returns
        -------
            val_out : array-like Quantity

        """

        if colname in self.meta:
            val = self.meta[colname]
        elif colname in self.table.colnames:
            val = self.table[colname]
        else:
            warnings.warn("{} not found in either '.meta' or '.table'"
                          "".format(colname))
            return None

        col_units = colname+"_unit"
        if isinstance(val, u.Quantity):
            units = val.unit
        elif col_units in self.meta:
            units = u.Unit(self.meta[col_units])
        else:
            units = u.Unit("")

        if isinstance(val, u.Quantity):
            val_out = val.to(units)
        elif isinstance(val, (list, tuple, np.ndarray)):
            val_out = val * units
        elif val is None:
            val_out = None
        else:
            raise ValueError("{} must be of type: Quantity, array, list, tuple"
                             "".format(colname))

        return val_out


def get_meta_quantity(meta_dict, name, fallback_unit=""):
    """
    Extract a Quantity from a dictionary

    Parameters
    ----------
    meta_dict : dict
    name : str
    fallback_unit : Quantity

    Returns
    -------
    quant : Quantity

    """

    if isinstance(meta_dict[name], u.Quantity):
        unit = meta_dict[name].unit
    elif name + "_unit" in meta_dict:
        unit = meta_dict[name + "_unit"]
    else:
        unit = u.Unit(fallback_unit)
    quant = quantify(meta_dict[name], unit)

    return quant


def quantify(item, unit):
    """
    Ensure an item is a Quantity

    Parameters
    ----------
    item : int, float, array, list, Quantity
    unit : str, Unit

    Returns
    -------
    quant : Quantity

    """

    if isinstance(item, u.Quantity):
        quant = item.to(u.Unit(unit))
    else:
        quant = item * u.Unit(unit)
    return quant


def make_emission_from_emissivity(temp, emiss_src_spec):
    """
    Create an emission SourceSpectrum using a blackbody and an emissivity curve

    Parameters
    ----------
    temp : float, Quantity
        [deg_C] If float, then must be in degrees Celsius
    emiss_src_spec : synphot.SpectralElement
        An emissivity response curve in the range [0..1]

    Returns
    -------
    flux : synphot.SourceSpectrum

    """

    if isinstance(temp, u.Quantity):
        temp = temp.to(u.deg_C)

    if emiss_src_spec is None:
        warnings.warn("Either emission or emissivity must be set")
        flux = None
    else:
        flux = SourceSpectrum(BlackBody1D, temperature=temp)
        flux.meta["solid_angle"] = u.sr**-1
        flux = flux * emiss_src_spec
        flux.meta["history"] = ["Created from Blackbody curve. Units are to be"
                                "understood as per steradian"]

    return flux


def make_emission_from_array(flux, wave, meta):
    """
    Create an emission SourceSpectrum using array.

    Takes care of bins and solid angles. The solid_angle is kept in the returned
    SourceSpectrum meta dictionary under self.meta["solid_angle"]

    Parameters
    ----------
    flux : array-like, Quantity
        if flux is not an array, the ``emission_unit`` must be in meta dict
    wave : array-like, Quantity
        if flux is not an array, the ``wavelength_unit`` must be in meta dict
    meta : dict

    Returns
    -------
    flux : synphot.SourceSpectrum

    """

    if not isinstance(flux, u.Quantity):
        if "emission_unit" in meta:
            flux = quantify(flux, meta["emission_unit"])
        else:
            warnings.warn("emission_unit must be set in self.meta, "
                          "or emission must be an astropy.Quantity")
            flux = None

    if isinstance(wave, u.Quantity) and isinstance(flux, u.Quantity):
        flux_unit, angle = extract_type_from_unit(flux.unit, "solid angle")
        flux = flux / angle

        if is_flux_binned(flux.unit):
            flux = normalise_binned_flux(flux, wave)

        orig_unit = flux.unit
        flux = SourceSpectrum(Empirical1D, points=wave,
                              lookup_table=flux)
        flux.meta["solid_angle"] = angle
        flux.meta["history"] = ["Created from emission array with units {}"
                                "".format(orig_unit)]
    else:
        warnings.warn("wavelength and emission must be "
                      "astropy.Quantity objects")
        flux = None

    return flux


def extract_type_from_unit(unit, unit_type):
    """
    Extract ``astropy`` physical type from a compound unit

    Parameters
    ----------
    unit : astropy.Unit
    unit_type : str
        The physical type of the unit as given by ``astropy``

    Returns
    -------
    new_unit : Unit
        The input unit minus any base units corresponding to ``unit_type``
    extracted_units : Unit
        Any base units corresponding to ``unit_type``

    """

    unit = unit**1
    extracted_units = u.Unit("")
    for base, power in zip(unit._bases, unit._powers):
        if unit_type == (base**abs(power)).physical_type:
            extracted_units *= base**power

    new_unit = unit / extracted_units

    return new_unit, extracted_units


def extract_base_from_unit(unit, base_unit):
    """
    Extract ``astropy`` base unit from a compound unit

    Parameters
    ----------
    unit : astropy.Unit
    base_unit : Unit, str

   Returns
    -------
    new_unit : Unit
        The input unit minus any base units corresponding to ``base_unit``
    extracted_units : Unit
        Any base units corresponding to ``base_unit``

    """

    unit = unit**1
    extracted_units = u.Unit("")
    for base, power in zip(unit._bases, unit._powers):
        if base == base_unit:
            extracted_units *= base**power

    new_unit = unit * extracted_units**-1

    return new_unit, extracted_units


def normalise_binned_flux(flux, wave):
    """
    Convert a binned flux Quantity array back into flux density

    The flux density normalising unit is taken from the wavelength Quantity unit

    Parameters
    ----------
    flux : array-like Quantity
        flux unit must include ``bin``, e.g. ``ph s-1 m-2 bin-1``
    wave : array-like Quantity

    Returns
    -------
    flux : array-like Quantity

    """

    bins = np.zeros(len(wave)) * wave.unit
    bins[:-1] = 0.5 * np.diff(wave)
    bins[1:] += 0.5 * np.diff(wave)
    # bins[0] *= 2.   # edge bins only have half the flux of other bins
    # bins[-1] *= 2.

    bin_unit = extract_base_from_unit(flux.unit, u.bin)[1]
    flux = flux / bins / bin_unit

    return flux


def is_flux_binned(unit):
    """
    Checks if the (flux) unit is a binned unit

    Parameters
    ----------
    unit : Unit

    Returns
    -------
    flag : bool

    """
    unit = unit**1
    flag = False
    if u.bin in unit._bases or "flux density" not in unit.physical_type:
        flag = True

    return flag
