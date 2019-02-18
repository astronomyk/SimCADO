import warnings
from copy import deepcopy

import numpy as np

from astropy import wcs as apwcs
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, Column

from . import image_plane_utils as imp_utils
from ..source.source2 import Source
from ..source import source2_utils as src_utils 

from .. import utils
from .. import rc


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

    def extract_from(self, src):
        """ ..assumption: Bandpass has been applied"""
        
        if not isinstance(src, Source):
            raise ValueError("source must be a Source object: {}"
                             "".format(type(src)))

        wave_min = utils.quantify(self.meta["wave_min"], u.um).value
        wave_max = utils.quantify(self.meta["wave_max"], u.um).value
        area = self.meta["area"]

        # determine which fields are inside the field of view
        fields_mask = [is_field_in_fov(self.hdu.header, field)
                       for field in src.fields]
        fields_indexes = np.where(fields_mask)[0]

        # make tables with x,y,flux for the relevant fields
        # determine which table rows are in FOV
        # get the photons_in_range for the relevant fields
        combined_table = combine_table_fields(self.hdu.header, src,
                                              fields_indexes)
        tbl = make_flux_table(combined_table, src, wave_max, wave_min, area)
        self.fields += [tbl]

        # make new imagehdu for the relevant fields
        # get photons_in_range for each image_spectrum
        # multiply the images by the flux
        # add_imagehdu_to_imagehdu

        imagehdu = combine_imagehdu_fields(self.hdu.header, src, fields_indexes,
                                           wave_max, wave_min, area)
        self.fields += [imagehdu]

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


def make_flux_table(source_tbl, src, wave_max, wave_min, area):
    fluxes = np.zeros(len(src.spectra))
    ref_set = list(set(source_tbl["ref"]))
    flux_set = src.photons_in_range(wave_min, wave_max, area, ref_set)
    fluxes[ref_set] = flux_set

    ref = source_tbl["ref"]
    weight = source_tbl["weight"]
    flux_col = Column(name="flux", data=fluxes[ref] * weight)
    x_col = source_tbl["x"]
    y_col = source_tbl["y"]

    tbl = Table()
    tbl.add_columns([x_col, y_col, flux_col])

    return tbl


def combine_table_fields(fov_header, src, field_indexes):
    fov_xsky, fov_ysky = imp_utils.calc_footprint(fov_header)

    x, y, ref, weight = [], [], [], []

    for ii in field_indexes:
        field = src.fields[ii]
        if isinstance(field, Table):
            xcol = utils.quantity_from_table("x", field, u.arcsec)
            ycol = utils.quantity_from_table("y", field, u.arcsec)
            x += list(xcol.to(u.deg).value)
            y += list(ycol.to(u.deg).value)
            ref += list(field["ref"])
            weight += list(field["weight"])

    x = np.array(x)
    y = np.array(y)
    mask = np.array(x < max(fov_xsky)) * np.array(x > min(fov_xsky)) * \
           np.array(y < max(fov_ysky)) * np.array(y > min(fov_ysky))

    x = x[mask]
    y = y[mask]
    ref = np.array(ref)[mask]
    weight = np.array(weight)[mask]

    tbl = Table(names=["x", "y", "ref", "weight"], data=[x, y, ref, weight])
    tbl["x"].unit = u.deg
    tbl["y"].unit = u.deg

    return tbl


def combine_imagehdu_fields(fov_header, src, fields_indexes,
                            wave_max, wave_min, area):
    image = np.zeros((fov_header["NAXIS1"], fov_header["NAXIS2"]))
    canvas_hdu = fits.ImageHDU(header=fov_header, data=image)
    order = int(rc.__rc__["SIM_SPLINE_ORDER"])

    for ii in fields_indexes:
        if isinstance(src.fields[ii], fits.ImageHDU):
            flux = src.photons_in_range(wave_min, wave_max, area, indexes=[ii])
            image = np.zeros((fov_header["NAXIS1"], fov_header["NAXIS2"]))
            temp_hdu = fits.ImageHDU(header=fov_header, data=image)
            temp_hdu = imp_utils.add_imagehdu_to_imagehdu(src.fields[ii],
                                                          temp_hdu, order=order)
            canvas_hdu.data += temp_hdu.data * flux[0].value

    return canvas_hdu

