import numpy as np
from scipy.signal import convolve

from astropy.io import fits

from ..fov import FieldOfView
from .effects import Effect


class PSF(Effect):
    def __init__(self, **kwargs):
        self.kernel = None
        self.valid_waverange = None
        super(Effect, self).__init__(**kwargs)

    def apply_to(self, fov):

        if not isinstance(fov, FieldOfView):
            raise ValueError("fov must be a FieldOfView object: {}"
                             "".format(type(fov)))
        if fov.hdu.data is None:
            fov.view()

        kernel = self.get_kernel(fov)

        old_shape = fov.hdu.data.shape
        new_image = convolve(fov.hdu.data, kernel, mode="full")
        new_shape = new_image.shape

        fov.hdu.data = new_image

        # ..todo: careful with which dimensions mean what
        fov.hdu.header["CRPIX1"] += (new_shape[0] - old_shape[0]) / 2
        fov.hdu.header["CRPIX2"] += (new_shape[1] - old_shape[1]) / 2
        fov.hdu.header["CRPIX1D"] += (new_shape[0] - old_shape[0]) / 2
        fov.hdu.header["CRPIX2D"] += (new_shape[1] - old_shape[1]) / 2

        return fov

    def fov_grid(self, header=None, waverange=None, **kwargs):
        return {"wavelengths": waverange}

    def get_kernel(self, fov):
        self.valid_waverange = None
        self.kernel = np.ones((1, 1))
        return self.kernel


################################################################################
# Analytical PSFs - Vibration, Seeing, NCPAs


class AnalyticalPSF(PSF):
    def __init__(self, **kwargs):
        super(PSF, self).__init__(**kwargs)


class Vibration(PSF):
    def __init__(self, **kwargs):
        super(PSF, self).__init__(**kwargs)

    def fov_grid(self, header=None, waverange=None, **kwargs):
        return {"coords": None, "wavelengths": None}


class NonCommonPathAberration(PSF):
    def __init__(self, **kwargs):
        super(PSF, self).__init__(**kwargs)

    def fov_grid(self, header=None, waverange=None, **kwargs):
        waveset = []
        return {"coords": None, "wavelengths": waveset}


class Seeing(PSF):
    def __init__(self, **kwargs):
        super(PSF, self).__init__(**kwargs)

    def fov_grid(self, header=None, waverange=None, **kwargs):
        waveset = []
        return {"coords": None, "wavelengths": waveset}


################################################################################
# Semi-analytical PSFs - Poppy PSFs


class SemiAnalyticalPSF(PSF):
    def __init__(self, **kwargs):
        super(PSF, self).__init__(**kwargs)


class PoppyFieldVaryingPSF(SemiAnalyticalPSF):
    def __init__(self, **kwargs):
        super(SemiAnalyticalPSF, self).__init__(**kwargs)


class PoppyFieldConstantPSF(SemiAnalyticalPSF):
    def __init__(self, **kwargs):
        super(SemiAnalyticalPSF, self).__init__(**kwargs)


################################################################################
# Discreet PSFs - MAORY and co PSFs


class DiscreetPSF(PSF):
    def __init__(self, **kwargs):
        super(PSF, self).__init__(**kwargs)


class MaoryFieldVaryingPSF(DiscreetPSF):
    def __init__(self, **kwargs):
        super(DiscreetPSF, self).__init__(**kwargs)


class MaoryFieldConstantPSF(DiscreetPSF):
    def __init__(self, **kwargs):
        super(DiscreetPSF, self).__init__(**kwargs)
        self.waveset, self.kernel_indexes = get_psf_wave_exts(self)
        self.current_layer = None

    def fov_grid(self, header=None, waverange=None, **kwargs):
        return {"wavelengths": self.waveset}

    def get_kernel(self, fov):
        fov_wave = 0.5 * (fov.meta["wave_min"] + fov.meta["wave_max"])
        ii = nearest_index(fov_wave, self.waveset)
        ext = self.kernel_indexes[ii]
        if ext != self.current_layer:
            self.kernel = self._file[ext]
            self.current_layer = ext


def nearest_index(x, x_array):
    return int(round(np.interp(x, x_array, np.arange(len(x_array)))))


def get_psf_wave_exts(hdu_list):
    """
    Returns a dict of {extension : wavelength}

    Parameters
    ----------
    hdu_list

    Returns
    -------
    wave_set, wave_ext

    """

    if not isinstance(hdu_list, fits.HDUList):
        raise ValueError("psf_effect must be a PSF object: {}"
                         "".format(type(hdu_list)))

    wave_set = [hdu.header["WAVE0"] for hdu in hdu_list
                if "WAVE0" in hdu.header]
    wave_ext = [ii for ii in range(len(hdu_list))
                if "WAVE0" in hdu_list[ii].header]

    return wave_set, wave_ext


