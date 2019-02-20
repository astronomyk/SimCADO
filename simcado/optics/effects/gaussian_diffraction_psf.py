import numpy as np
from astropy import units as u
from astropy.convolution import Gaussian2DKernel
from scipy.signal import fftconvolve

from ... import utils
from ..fov import FieldOfView
from .effects import Effect


class GaussianDiffractionPSF(Effect):
    def __init__(self, diameter, **kwargs):
        super(GaussianDiffractionPSF, self).__init__(**kwargs)
        self.meta["diameter"] = diameter

    def apply_to(self, fov, **kwargs):
        if not isinstance(fov, FieldOfView):
            raise ValueError("fov must be a FieldOfView object: {}"
                             "".format(type(fov)))

        pixel_scale = fov.header["CDELT1"] * u.deg.to(u.arcsec)
        wave = 0.5 * (fov.meta["wave_max"] + fov.meta["wave_min"])
        kernel = get_diffraction_limited_gaussian_kernel(wave, pixel_scale,
                                                         self.meta["diameter"])

        new_im = fftconvolve(fov.data, kernel, mode="full")
        fov.hdu.data = new_im
        fov.hdu.header["CRPIX1"] += kernel.shape[0] / 2.
        fov.hdu.header["CRPIX2"] += kernel.shape[0] / 2.

        return fov

    def fov_grid(self, header, waverange):
        waverange = utils.quantify(waverange, u.um)
        diameter = utils.quantify(self.meta["diameter"], u.m).to(u.um)
        fwhm = 1.22 * (waverange / diameter).value  # in rad

        pixel_scale = utils.quantify(header["CDELT1"], u.deg).to(u.rad).value
        fwhm_range = np.arange(fwhm[0], fwhm[1], pixel_scale)
        wavelengths = fwhm_range / 1.22 * diameter.to(u.m)

        return {"coords": None, "wavelengths": wavelengths}

    def update(self, **kwargs):
        if "diameter" in kwargs:
            self.meta["diameter"] = kwargs["diameter"]


def get_diffraction_limited_gaussian_kernel(wave, pixel_scale, diameter):
    wave = utils.quantify(wave, u.um)
    diameter = utils.quantify(diameter, u.m).to(u.um)
    pixel_scale = utils.quantify(pixel_scale, u.arcsec)
    fwhm = 1.22 * (wave / diameter) * u.rad.to(u.arcsec) / pixel_scale

    sigma = fwhm.value / 2.35
    kernel = Gaussian2DKernel(sigma, mode="center").array

    return kernel
