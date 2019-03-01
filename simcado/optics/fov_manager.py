import numpy as np
from astropy import units as u

from . import effects as efs
from .fov import FieldOfView
from .image_plane_utils import header_from_list_of_xy
from .effects.effects_utils import _get_all, is_spectroscope


# 1. Find the Wavelength range
# Build from edges of throughput curve

# 2. Find the wavelength bins
# If TraceList and Aperture list, then Spectroscopy
# TraceList
# for each trace dlam along the trace centre in increments
#   of SIM_SUB_PIXEL_FRACTION
# Must be accompanied by an ApertureList

# If not, then imaging
# PSF core increase (atmo, ncpas)
# If from a files, what is the bin size?
# If analytic, dlam between a FWHM or SIM_SUB_PIXEL_FRACTION
# ADC + AD shifts
# dlam between shift of SIM_SUB_PIXEL_FRACTION

# 3. Find the spatial range
# If Spectroscopy
# ApertureList
# For each Trace set the sky header to the aperture footprint
#   plus any shifts from AtmosphericDispersion
# Set the Image plane footprint centred on the image plane
#   position

# If Imaging
# DetectorList, or ApertureMask, plus any shift from
#   AtmosphericDispersion


class FOVManager:
    def __init__(self, effects=[], **kwargs):
        self.meta = {}
        self.meta.update(kwargs)
        self.effects = effects
        self._fovs_list = []

    def generate_fovs_list(self):
        waverange = [self.meta["SIM_LAM_MIN"] * u.um,
                     self.meta["SIM_LAM_MAX"] * u.um]
        waverange = self.effects[0].fov_grid(None, waverange)["wavelength"]

        if is_spectroscope(self.effects):
            shifts  = get_3D_shifts(self.effects, **self.meta)
            waveset = get_spectroscopy_waveset(self.effects, **self.meta)
            headers = get_spectroscopy_headers(self.effects, **self.meta)
            hdus    = get_spectroscopy_hdus(headers, waveset, shifts)
        else:
            shifts  = get_3D_shifts(self.effects, **self.meta)
            waveset = get_imaging_waveset(self.effects, **self.meta)
            headers = get_imaging_headers(self.effects, **self.meta)
            hdus    = get_imaging_hdus(headers, waveset, shifts)

        return hdus

    @property
    def fovs(self):
        self._fovs_list = self.generate_fovs_list()
        return self._fovs_list


def get_3D_shifts(effects, **kwargs):
    # for the offsets
    # OBS_FIELD_ROTATION
    # ATMO_TEMPERATURE
    # ATMO_PRESSURE
    # ATMO_REL_HUMIDITY
    # ATMO_PWV
    # ATMO_AIRMASS

    wave_min = kwargs["SIM_LAM_MIN"]
    wave_mid = kwargs["SIM_LAM_MID"]
    wave_max = kwargs["SIM_LAM_MAX"]
    sub_pixel_frac = kwargs["SIM_SUB_PIXEL_FRACTION"]

    effects = _get_all(effects, efs.Shift3D)
    shifts = [eff.fov_grid(lam_min=wave_min, lam_mid=wave_mid, lam_max=wave_max,
                           sub_pixel_frac=sub_pixel_frac, **kwargs)
              for eff in effects]

    # ..todo: Set this up so that it actually does something
    wave_bin_edges = [wave_min, wave_max]
    x_shifts = [0, 0]
    y_shifts = [0, 0]

    return {"wavelengths": wave_bin_edges,
            "x_shifts": x_shifts,
            "y_shifts": y_shifts}


def get_imaging_waveset(effects, **kwargs):
    # dlam between shift of

    wave_min = kwargs["SIM_LAM_MIN"]
    wave_mid = kwargs["SIM_LAM_MID"]
    wave_max = kwargs["SIM_LAM_MAX"]
    sub_pixel_frac = kwargs["SIM_SUB_PIXEL_FRACTION"]

    psfs = _get_all(effects, efs.PSF)
    if len(psfs) > 0:
        wave_bin_edges = [psf.fov_grid(sub_pixel_frac=sub_pixel_frac,
                                       lam_min=wave_min,
                                       lam_max=wave_max)["wavelengths"]
                          for psf in psfs]
        # assume the longest array requires the highest spectral resolution
        len_steps = [len(lbe) for lbe in wave_bin_edges]
        ii = np.where(len_steps == max(len_steps))[0][0]
        wave_bin_edges = wave_bin_edges[ii]
    else:
        wave_bin_edges = [wave_min, wave_max]

    return wave_bin_edges


def get_imaging_headers(effects, waveset, **kwargs):

    # check DetectorList for boundaries, convert to on-sky coords
    # check ApertureMask for boundaries
    # find the smallest region
    # cut up based on MAX_CHUNK_SIZE

    aperture_masks = _get_all(effects, efs.ApertureMask)
    detector_array = _get_all(effects, efs.DetectorList)[0]

    pixel_size = detector_array.image_plane_header["CDELT1"]
    pixel_scale = kwargs["SIM_DETECTOR_PIX_SCALE"]
    deg2mm = pixel_size / pixel_scale

    if len(aperture_masks) > 0:
        sky_edges = [apm.fov_grid()["edges"] for apm in aperture_masks]
    elif len(detector_array) > 0:
        sky_edges = [detector_array.fov_grid(pixel_scale=pixel_scale)["edges"]]
    else:
        raise ValueError("No ApertureMask or DetectorList was provided. At "
                         "least a DetectorList object must be passed: {}"
                         "".format(effects))

    width = kwargs["SIM_CHUNK_SIZE"]
    hdrs = []
    for xy_sky in sky_edges:
        x0, y0 = min(xy_sky[0]), min(xy_sky[1])
        x1, y1 = max(xy_sky[0]), max(xy_sky[1])
        for xi in np.arange(x0, x1, width):
            for yi in np.arange(y0, y1, width):
                xii = np.array([xi, xi + min(width, x1-xi)])
                yii = np.array([yi, yi + min(width, y1-yi)])
                hdr_sky = header_from_list_of_xy(xii, yii, pixel_scale)
                hdr_mm  = header_from_list_of_xy(xii * deg2mm, yii * deg2mm,
                                                 pixel_size, "D")
                hdr_sky.update(hdr_mm)
                hdrs += [hdr_sky]

    return hdrs


def get_imaging_hdus(headers, waveset, shifts):

    if len(shifts["wavelength"]) > len(waveset):
        waveset = shifts["wavelength"]

    # ..todo: add the shifts in somehow

    fovs = []
    for ii in range(len(waveset) - 1):
        for hdr in headers:
            waverange = [waveset[ii], waveset[ii + 1]]
            fovs += [FieldOfView(hdr, waverange)]

    return fovs


def get_spectroscopy_waveset(effects, **kwargs):
    # TraceList
    # for each trace dlam along the trace centre in increments
    #   of SIM_SUB_PIXEL_FRACTION
    # Must be accompanied by an ApertureList
    pass


def get_spectroscopy_headers(effects, **kwargs):
    # ApertureList
    # For each Trace set the sky header to the aperture footprint
    #   plus any shifts from AtmosphericDispersion
    # Set the Image plane footprint centred on the image plane
    #   position
    pass


def get_spectroscopy_hdus(fields, waveset, shifts):
    pass