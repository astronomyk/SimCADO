import warnings
import numpy as np

from astropy import wcs as apwcs
from astropy import units as u
from astropy.io import fits
from astropy.table import Table
from astropy import wcs

from synphot.units import PHOTLAM

from . import image_plane_utils as imp_utils
from ..source.source2 import Source
from ..source import source2_utils as src_utils 

from .. import utils
from .. import rc


# WARNING : all spectra should be resampled so that there are at least 10 bins
# per wavelength layer before we start for the FOV-ing


class FieldOfView:
    """
    A FOV is a monochromatic image. Flux units after extracting the fields from
    the Source are in ph/s/pixel

    The initial header should contain a spatial WCS including:
    - CDELTn, CUNITn, NAXISn : for pixel scale and size (assumed CUNIT in deg)
    - CRVALn, CRPIXn : for positioning the final image
    - CTYPE : is assumed to be a TAN projection.


    """

    def __init__(self, header, waverange, **kwargs):
        self.meta = {"wave_binwidth" : rc.__rc__["SIM_SPEC_RESOLUTION"],
                     "wave_min" : utils.quantify(waverange[0], u.um),
                     "wave_max" : utils.quantify(waverange[1], u.um),
                     "area" : 1 * u.m**2,
                     "sub_pixel" : rc.__rc__["SIM_SUB_PIXEL_ACCURACY"]}
        self.meta.update(kwargs)

        if len(apwcs.find_all_wcs(header)) == 0:
            raise ValueError("header must contain a valid WCS: {}"
                             "".format(dict(header)))

        data = np.zeros((header["NAXIS1"], header["NAXIS2"]))
        self.hdu = fits.ImageHDU(header=header, data=data)
        self.fields = []
        self.fluxes = []

    def extract_from(self, source):
        """ ..assumption: Bandpass has been applied"""
        
        if not isinstance(source, Source):
            raise ValueError("source must be a Source object: {}"
                             "".format(type(source)))

        # determine which fields are inside the field of view
        # determine which table rows are in FOV
        # get the photons_in_range for the relevant fields
        # make tables with x,y,flux for the relevant fields
        # make new imagehdu for the relevant fields

        fields_mask = [is_field_in_fov(self.hdu.header, field)
                       for field in source.fields]
        fields_indexes = np.where(fields_mask)[0]

        print(fields_indexes)

        wave_min = utils.quantify(self.meta["wave_min"], u.um).value
        wave_max = utils.quantify(self.meta["wave_max"], u.um).value
        #
        # self.fields += [make_flux_table(self.hdu.header, source, fields_indexes)]
        # self.fields += [reduce_imagehdus(self.fields, self.spectra,
        #                                     fields_indexes)]
        #
        #
        #








        # fluxes = []
        # int_flux_unit = PHOTLAM * u.Angstrom
        # for spec in source.spectra:
        #     wave = spec.model.points[0]
        #     flux = spec.model.lookup_table
        #
        #     mask = (spec.model.points[0] >= wave_min) * \
        #            (spec.model.points[0] <= wave_max)
        #     fluxes += [np.trapz(flux[mask], wave[mask])]
        #
        # fluxes = fluxes * int_flux_unit
        #
        # for field in source.fields:
        #     if isinstance(field, Table):
        #         x = utils.quantity_from_table("x", field, u.arcsec)
        #         y = utils.quantity_from_table("y", field, u.arcsec)
        #         flux = fluxes[tbl["ref"]]
        #         tbl = Table(names=["x", "y", "flux"], data=[x, y, flux])
        #
        #         sub_pixel = self.meta["sub_pixel"]
        #         self.hdu = imp_utils.add_table_to_imagehdu(tbl, self.hdu,
        #                                                     sub_pixel=sub_pixel)
        #     elif isinstance(field, fits.ImageHDU):
        #         if field.header["SPEC_REF"] is not "":
        #             flux =  fluxes[field.header["SPEC_REF"]]
        #         else:
        #             flux = 1 * int_flux_unit
        #
        #         self.hdu = imp_utils.add_imagehdu_to_imagehdu(field, self.hdu)


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


def is_field_in_fov(fov_header, table_or_imagehdu):

    pixel_scale = utils.quantify(fov_header["CDELT1"], u.deg)

    if isinstance(table_or_imagehdu, Table):
        ext_hdr = imp_utils._make_bounding_header_for_tables(
                                            [table_or_imagehdu], pixel_scale)
    elif isinstance(table_or_imagehdu, fits.ImageHDU):
        ext_hdr = imp_utils._make_bounding_header_from_imagehdus(
                                            [table_or_imagehdu], pixel_scale)
    else:
        raise ValueError(table_or_imagehdu)

    ext_xsky, ext_ysky = imp_utils.calc_footprint(ext_hdr)
    fov_xsky, fov_ysky = imp_utils.calc_footprint(fov_header)

    is_inside_fov = min(ext_xsky) < max(fov_xsky) and \
                    max(ext_xsky) > min(fov_xsky) and \
                    min(ext_ysky) < max(fov_ysky) and \
                    max(ext_ysky) > min(fov_ysky)

    return is_inside_fov


def make_flux_table(fov_header, src, indexes):
    fov_xsky, fov_ysky = imp_utils.calc_footprint(fov_header)

    x, y, ref, weight = [], [], [], []
    for ii, field in enumerate(src.fields):
        if isinstance(field, Table) and indexes[ii] is True:
            xcol = utils.quantity_from_table(field, "x", u.arcsec)
            ycol = utils.quantity_from_table(field, "y", u.arcsec)
            x += list(xcol.to(u.deg).value)
            y += list(ycol.to(u.deg).value)
            ref += list(field["ref"])
            weight += list(field["weight"])

    x = np.array(x)
    y = np.array(y)
    mask = x < max(fov_xsky) * x > min(fov_xsky) * \
           y < max(fov_ysky) * y > min(fov_ysky)

    spec_refs = set(ref[mask])