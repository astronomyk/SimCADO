import warnings
import numpy as np

from astropy import wcs as apwcs
from astropy import units as u
from astropy.io import fits
from astropy.table import Table

from synphot.units import PHOTLAM

from . import image_plane_utils as impl_utils

from ..source.source2 import Source
from ..import utils
from .. import rc


# WARNING : all spectra should be resampled so that there are at least 10 bins
# per wavelength layer before we start for the FOV-ing


class FieldOfView:
    def __init__(self, header, waverange, **kwargs):
        self.meta = {"wave_binwidth" : rc.__rc__["SIM_SPEC_RESOLUTION"],
                     "wave_min" : utils.quantify(waverange[0], u.um),
                     "wave_max" : utils.quantify(waverange[1], u.um),
                     "area" : 1 * u.m**2,
                     "transmission" : None,
                     "sub_pixel" : rc.__rc__["SIM_SUB_PIXEL_ACCURACY"]}
        self.meta.update(kwargs)

        self.hdu = fits.ImageHDU(header=header)
        if len(apwcs.find_all_wcs(header)) == 0:
            raise ValueError("header must contain a valid WCS: {}"
                             "".format(dict(header)))

    def extract_from(self, source):
        if not isinstance(source, Source):
            raise ValueError("source must be a Source object: {}"
                             "".format(type(source)))

        wave_min = self.meta["wave_min"].to(u.Angstrom).value
        wave_max = self.meta["wave_max"].to(u.Angstrom).value
        fluxes = []
        int_flux_unit = PHOTLAM * u.Angstrom
        for spec in source.spectra:
            wave = spec.model.points[0]
            flux = spec.model.lookup_table

            mask = (spec.model.points[0] >= wave_min) * \
                   (spec.model.points[0] <= wave_max)
            fluxes += [np.trapz(flux[mask], wave[mask])]

        fluxes = fluxes * int_flux_unit

        for field in source.fields:
            if isinstance(field, Table):
                x = utils.quantity_from_table("x", field, u.arcsec)
                y = utils.quantity_from_table("y", field, u.arcsec)
                flux = fluxes[tbl["ref"]]
                tbl = Table(names=["x", "y", "flux"], data=[x, y, flux])

                sub_pixel = self.meta["sub_pixel"]
                self.hdu = impl_utils.add_table_to_imagehdu(tbl, self.hdu,
                                                            sub_pixel=sub_pixel)
            elif isinstance(field, fits.ImageHDU):
                if field.header["SPEC_REF"] is not "":
                    flux =  fluxes[field.header["SPEC_REF"]]
                else:
                    flux = 1 * int_flux_unit

                self.hdu = impl_utils.add_imagehdu_to_imagehdu(field, self.hdu)


    ###################################################
    # Deal with flux units properly!!!!!

    @property
    def header(self):
        return self.hdu.header

    @property
    def data(self):
        return self.hdu.data

    @property
    def image(self):
        return self.data
