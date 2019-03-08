import numpy as np
from scipy.signal import convolve

from astropy.io import fits

from ..fov import FieldOfView
from .effects import Effect


class PSF(Effect):
    def __init__(self, **kwargs):
        self.kernel = None
        self.valid_waverange = None
        super(PSF, self).__init__(**kwargs)

    def apply_to(self, obj):
        if len(obj.fields) > 0:
            kernel = self.get_kernel(obj)

            sub_pixel = False
            if "SIM_SUB_PIXEL_FLAG" in self.meta:
                sub_pixel = self.meta["SIM_SUB_PIXEL_FLAG"]

            if obj.hdu.data is None:
                obj.view(sub_pixel)

            old_shape = obj.hdu.data.shape
            new_image = convolve(obj.hdu.data, kernel, mode="full")
            new_shape = new_image.shape

            obj.hdu.data = new_image

            # ..todo: careful with which dimensions mean what
            obj.hdu.header["CRPIX1"] += (new_shape[0] - old_shape[0]) / 2
            obj.hdu.header["CRPIX2"] += (new_shape[1] - old_shape[1]) / 2
            obj.hdu.header["CRPIX1D"] += (new_shape[0] - old_shape[0]) / 2
            obj.hdu.header["CRPIX2D"] += (new_shape[1] - old_shape[1]) / 2

        return obj

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
        super(AnalyticalPSF, self).__init__(**kwargs)


class Vibration(AnalyticalPSF):
    def __init__(self, **kwargs):
        super(Vibration, self).__init__(**kwargs)


class NonCommonPathAberration(AnalyticalPSF):
    def __init__(self, **kwargs):
        super(NonCommonPathAberration, self).__init__(**kwargs)


class Seeing(AnalyticalPSF):
    def __init__(self, **kwargs):
        super(Seeing, self).__init__(**kwargs)


################################################################################
# Semi-analytical PSFs - Poppy PSFs


class SemiAnalyticalPSF(PSF):
    def __init__(self, **kwargs):
        super(SemiAnalyticalPSF, self).__init__(**kwargs)


class PoppyFieldVaryingPSF(SemiAnalyticalPSF):
    def __init__(self, **kwargs):
        super(PoppyFieldVaryingPSF, self).__init__(**kwargs)


class PoppyFieldConstantPSF(SemiAnalyticalPSF):
    def __init__(self, **kwargs):
        super(PoppyFieldConstantPSF, self).__init__(**kwargs)


################################################################################
# Discreet PSFs - MAORY and co PSFs


class DiscretePSF(PSF):
    def __init__(self, **kwargs):
        super(DiscretePSF, self).__init__(**kwargs)


class MaoryFieldVaryingPSF(DiscretePSF):
    def __init__(self, **kwargs):
        super(MaoryFieldVaryingPSF, self).__init__(**kwargs)


class MaoryFieldConstantPSF(DiscretePSF):
    def __init__(self, **kwargs):
        super(MaoryFieldConstantPSF, self).__init__(**kwargs)
        self.waveset, self.kernel_indexes = get_psf_wave_exts(self)
        self.current_layer_id = None

    def fov_grid(self, header=None, waverange=None, **kwargs):
        return {"wavelengths": self.waveset}

    def get_kernel(self, fov):
        fov_wave = 0.5 * (fov.meta["wave_min"] + fov.meta["wave_max"])
        ii = nearest_index(fov_wave, self.waveset)
        ext = self.kernel_indexes[ii]
        if ext != self.current_layer_id:
            self.kernel = self._file[ext]
            self.current_layer_id = ext


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


