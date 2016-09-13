# pylint: disable=too-many-lines
"""
The module which contains the functionality to create Source objects

Summary


Classes
-------
Source

Functions
---------
Functions to create `Source` objects
```
empty_sky()
star(mag, filter_name="K", spec_type="A0V", position=(0,0))
stars(mags, x, y, filter_name="K", spec_types="A0V")
star_grid(n, mag_min, mag_max, filter_name="K",0 seperation=1, area=1,
          spec_type="A0V")
source_from_image(images, lam, spectra, pix_res, oversample=1, units="ph/s/m2",
                  flux_threshold=0, center_pixel_offset=(0,0))
source_1E4_Msun_cluster(distance=50000, half_light_radius=1)
```

Functions for manipulating spectra for a `Source` object
```
scale_spectrum(lam, spec, mag, filter_name="K", return_ec=False)
scale_spectrum_sb(lam, spec, mag_per_arcsec, pix_res=0.004, filter_name="K",
                      return_ec=False)
flat_spectrum(mag, filter_name="K", return_ec=False)
flat_spectrum_sb(mag_per_arcsec, filter_name="K", pix_res=0.004, return_ec=False)
```

Functions regarding photon flux and magnitudes
```
zero_magnitude_photon_flux(filter_name, area=1)
_get_stellar_properties(spec_type, cat=None, verbose=False)
_get_stellar_mass(spec_type)
_get_stellar_Mv(spec_type)
_get_pickles_curve(spec_type, cat=None, verbose=False)
```

Helper functions
```
value_at_lambda(lam_i, lam, val, return_index=False)
SED(spec_type, filter_name="V", magnitude=0.)
```

See also
--------

Examples
--------


"""
###############################################################################
# source
#
# DESCRIPTION
# The source is essentially a list of spectra and a list of positions. The
# list of positions contains a reference to the relevant spectra. The advantage
# here is that if there are repeated spectra in a data cube, we can reduce the
# amount of calculations needed. Furthermore, if the input is originally a list
# of stars, etc where the position of a star is not always and integer multiple
# of the plate scale, we can keep the information until the PSFs are needed.
#
# The source contains two arrays:
#  - PositionArray:
#  - SpectrumArray


# Flow of events
# - Generate the lists of spectra and positions
# - Apply the transmission curves [SpectrumArray]
# - shrink the 1D spectra to the resolution of the psf layers [SpectrumArray]
# - apply any 2D spatials [PositionArray]
# for i in len(slices)
#   - generate a working slice [PositionArray, SpectrumArray, WorkingSlice]
#   - Apply the PSF for the appropriate wavelength [WorkingSlice]
#   - Apply any wavelength dependent spatials [WorkingSlice]
#   - apply Poisson noise to the photons in the slice [WorkingSlice]
#   - add the WorkingSlice to the FPA [WorkingSlice, FPArray]


## TODO implement conversions to Source object from:
# ascii
#    x,y,mag,[temp]
#    x,y,type
# images
#    JHK
#    cube


import os
import inspect
__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))

from copy import deepcopy

import numpy as np
import scipy.ndimage.interpolation as spi

from astropy.io import fits
from astropy.io import ascii as ioascii
from astropy.convolution import convolve, convolve_fft
import astropy.units as u
import astropy.constants as c


try:
    import simcado.spatial as pe
    import simcado.spectral as sc
    import simcado.psf as sim_psf
    import simcado.utils as utils
except ImportError:
    import spatial as pe
    import spectral as sc
    import psf as sim_psf
    import utils

__pkg_dir__ = os.path.split(__file__)[0]
__all__ = ["Source"]


# add_uniform_background() moved to detector
# get_slice_photons() renamed to photons_in_range() and moved to Source
# apply_transmission_curve() moved to Source
# apply_optical_train() moved to Source


class Source(object):
    """
    Create a source object from a file or from arrays

    Source class generates the arrays needed for source. It takes various
    inputs and converts them to an array of positions and references to spectra
    It also converts spectra to photons/s/voxel. The default units for input
    data is ph/s

    Parameters
    ==========
    - filename
    or
    - lam
    - spectra
    - x [arcsec]
    - y [arcsec]
    - ref
    - weight

    Keyword arguments
    =================
    - units
    - pix_unit
    - exptime
    - area
    """

    def __init__(self, filename=None,
                 lam=None, spectra=None, x=None, y=None, ref=None, weight=None,
                 **kwargs):

        self.params = {"units"   :"ph/s",
                       "pix_unit":"arcsec",
                       "exptime" :1,
                       "area"    :1,
                       "pix_res" :0.004}
        self.params.update(kwargs)

        self.info = dict([])
        self.info['created'] = 'yes'
        self.info['description'] = "List of spectra and their positions"

        self.units = u.Unit(self.params["units"])
        self.exptime = self.params["exptime"]
        self.area = self.params["area"]
        self.pix_res = self.params["pix_res"]

        self.x = None
        self.y = None  # set later

        if filename is not None:
            hdr = fits.getheader(filename)
            if "SIMCADO" in hdr.keys() and hdr["SIMCADO"] == "SOURCE":
                self.read(filename)
            else:
                self._from_cube(self, filename)
        elif not None in (lam, spectra, x, y, ref):
            self._from_arrays(lam, spectra, x, y, ref, weight)
        else:
            raise ValueError("Trouble with inputs. Could not create Source")

        self.ref = np.array(self.ref, dtype=int)
        self.x_orig = deepcopy(self.x)
        self.y_orig = deepcopy(self.y)


    def apply_optical_train(self, opt_train, detector, chips="all", **kwargs):
        """
        Apply all effects along the optical path to the source photons

        Parameters
        ----------
        opt_train : simcado.OpticalTrain
            the object containing all information on what is to happen to the
            photons as they travel from the source to the detector
        detector : simcado.Detector
            the object representing the detector
        chips : int, list
            Default is "all"

        Notes
        -----
        Output array is in units of [ph/s/pixel] where the pixel is internal
        oversampled pixels - not the pixel size of the detector chips

        """
        params = {"verbose"     :opt_train.cmds.verbose,
                  "INST_DEROT_PERFORMANCE"  :100,
                  "SCOPE_JITTER_FWHM"       :0,
                  "SCOPE_DRIFT_DISTANCE"    :0}
        params.update(self.params)
        params.update(kwargs)

        self.pix_res = opt_train.pix_res

        # 1. Apply the master transmission curve to all the spectra
        #
        # 1.5 Create a canvas onto which we splat the PSFed sources
        #
        # 2. For each layer between cmds.lam_bin_edges[i, i+1]
        #   - Apply the x,y shift for the ADC
        #       - Apply any other shifts
        #   - Apply the PSF for the layer
        #       - Sum up the photons in this section of the spectra
        #   - Add the layer to the final image array
        #
        # 3. Apply wave-indep psfs
        #   - field rotation
        #   - telescope shake
        #   - tracking error
        #
        # 3.5 Up until now everything is ph/s/m2/bin
        #     Apply the collecting area of the telescope
        #
        # 4. Add the average number of atmo-bg and mirror-bb photons
        # 5. Apply the instrumental distortion

        if chips is None or str(chips).lower() == "all":
            chips = np.arange(len(detector.chips))

        if not hasattr(chips, "__len__"):
            chips = [chips]

        # 1.
        self.apply_transmission_curve(opt_train.tc_source)

        for chip_i in chips:

            # 1.5
            image = None

            # 2.
            for i in range(len(opt_train.lam_bin_edges[:-1])):

                if params["verbose"]:
                    print("Wavelength slice [um]:", \
                                            opt_train.lam_bin_centers[i])
                # apply the adc shifts
                self.x = self.x_orig + opt_train.adc_shifts[0][i]
                self.y = self.y_orig + opt_train.adc_shifts[1][i]

                # include any other shifts here


                # apply the psf (get_slice_photons is called within)
                lam_min, lam_max = opt_train.lam_bin_edges[i:i+2]
                psf = opt_train.psf[i]


                oversample = opt_train.cmds["SIM_OVERSAMPLING"]
                if image is None:
                    image = self.image_in_range(psf, lam_min, lam_max,
                                                detector.chips[chip_i],
                                                pix_res=opt_train.pix_res,
                                                oversample=oversample)
                else:
                    image += self.image_in_range(psf, lam_min, lam_max,
                                                 detector.chips[chip_i],
                                                 pix_res=opt_train.pix_res,
                                                 oversample=oversample)

            # 3.
            # !!!!!!!!!!!!!! All of these need to be combined into a single
            # function that traces out the path taken by the telescope,
            # rather than having the arcs from the derotator() function
            # being stretched by the tracking() function and then the whole
            # thing blurred by wind_jitter()
            if params["INST_DEROT_PERFORMANCE"] < 100:
                image = opt_train.apply_derotator(image)
            if params["SCOPE_DRIFT_DISTANCE"] > 0.1 * self.pix_res:
                image = opt_train.apply_tracking(image)
            if params["SCOPE_JITTER_FWHM"] > 0.1 * self.pix_res:
                image = opt_train.apply_wind_jitter(image)

            # 3.5
            image *= opt_train.cmds.area

            # 4.
            image += (opt_train.n_ph_atmo + opt_train.n_ph_mirror)
            ## TODO: protected members should not be set by another class (OC)
            detector._n_ph_atmo = opt_train.n_ph_atmo
            detector._n_ph_mirror = opt_train.n_ph_mirror
            # 5.
            self.project_onto_chip(image, detector.chips[chip_i])

        ######################################
        # CAUTION WITH THE PSF NORMALISATION #
        ######################################


    def project_onto_chip(self, image, chip):
        """
        Re-project the photons onto the same grid as the detectors use

        Parameters
        ----------
        image : np.ndarray
            the image to be re-projected
        chip : detector.Chip
            the chip object where the image will land
        """

        chip.reset()
        scale_factor = self.pix_res / chip.pix_res

        chip_arr = spi.zoom(image, scale_factor, order=1)
        chip_arr *= np.sum(image) / np.sum(chip_arr)

        chip.add_signal(chip_arr)


    def image_in_range(self, psf, lam_min, lam_max, chip, **kwargs):
        """
        Find the sources that fall in the chip area and generate an image for
        the wavelength range [lam_min, lam_max)

        Output is in [ph/s/pixel]

        Parameters
        ----------
        psf : psf.PSF object
            The PSF that the sources will be convolved with
        lam_min, lam_max : float
            [um] the wavelength range relevant for the psf
        chip : detector.Chip object
            the chip that will be seeing this image.

        Optional parameters (**kwargs)
        ------------------------------
        sub_pixel : bool
            if sub-pixel accuracy is needed, each source is shifted individually.
            Default is False
        pix_res : float
            [arcsec] the field of view of each pixel. Default is 0.004 arcsec
        oversample : int
            the psf images will be oversampled to better conserve flux.
            Default is 1 (i.e. not oversampled)

        Returns
        -------
        slice_array : np.ndarray
            the image of the source convolved with the PSF for the given range

        """

        params = {"pix_res"     :0.004,
                  "sub_pixel"   :False,
                  "oversample"  :1}

        params.update(kwargs)

        if isinstance(psf, (sim_psf.PSFCube, sim_psf.UserPSFCube)):
            lam_cen = (lam_max + lam_min) / 2.
            psf = psf.nearest(lam_cen)

        if isinstance(psf, np.ndarray):
            arr = deepcopy(psf)
            pix_res = params["pix_res"] / params["oversample"]
            size = psf.shape[0]
            psf = sim_psf.PSF(size, pix_res)
            psf.set_array(arr)


        if chip is not None:
            mask = (self.x > chip.x_min) * (self.x < chip.x_max) * \
                   (self.y > chip.y_min) * (self.y < chip.y_max)
            params["pix_res"] = chip.pix_res / params["oversample"]
            x_min, x_max = chip.x_min, chip.x_max,
            y_min, y_max = chip.y_min, chip.y_max
            x_cen, y_cen = chip.x_cen, chip.y_cen
        else:
            mask = np.array([True]*len(self.x))
            params["pix_res"] /= params["oversample"]
            x_min, x_max = np.min(self.x), np.max(self.x)
            y_min, y_max = np.min(self.y), np.max(self.y),
            x_cen, y_cen = (x_max + x_min) / 2, (y_max + y_min) / 2


        naxis1 = int((x_max - x_min) / params["pix_res"])
        naxis2 = int((y_max - y_min) / params["pix_res"])

        slice_array = np.zeros((naxis1, naxis2), dtype=np.float32)
        slice_photons = self.photons_in_range(lam_min, lam_max, min_bins=10)

        x_pix = (self.x - x_cen) / params["pix_res"]
        y_pix = (self.y - y_cen) / params["pix_res"]


        # if sub pixel accuracy is needed, be prepared to wait. For this we
        # need to go through every source spectra in turn, shift the psf by
        # the decimal amount given by pos - int(pos), then place a
        # certain slice of the psf on the output array.
        ax, ay = np.array(slice_array.shape) // 2
        bx, by = np.array(psf.array.shape)   // 2

        print("Creating layer between [um]:", lam_min, lam_max)

        if params["sub_pixel"] is True:
            # for each point source in the list, add a psf to the slice_array
            x_int, y_int = x_pix.astype(int), y_pix.astype(int)
            dx, dy = self.x - x_int, self.y - y_int

            for i in range(len(slice_photons)):

                if not mask[i]:
                    continue

                ## Unused variable (OC)
                # psf_tmp = spi.shift(psf.array, (dx[i], dy[i]), order=1)
                x_pint, y_pint = x_int[i], y_int[i]

                # Find the slice borders for the array where the psf will go
                ax0 = np.max(np.array((x_pint - bx, [0]*len(x_pint))), axis=0)
                ax1 = np.min(np.array((x_pint + bx + 1,
                                       [slice_array.shape[0]]*len(x_pint))),
                             axis=0)
                ay0 = np.max(np.array((y_pint - by, [0] * len(y_pint))),
                             axis=0)
                ay1 = np.min(np.array((y_pint + by + 1,
                                       [slice_array.shape[1]]*len(y_pint))),
                             axis=0)

                # the slice limits of the psf array are found by taking the
                # pixel distance from the x,y position to the slice limits
                # of the slice_array. This distance is subtracted from the
                # centre of the psf array.
                bx0 = bx - (x_pint - ax0)
                bx1 = bx + (ax1 - x_pint)
                by0 = by - (y_pint - ay0)
                by1 = by + (ay1 - y_pint)

                slice_array[ax0:ax1, ay0:ay1] = psf.array[bx0:bx1, by0:by1] \
                                        * slice_photons[i] * self.weight[i]

        elif params["sub_pixel"] == "raw":
            #x_int, y_int = np.round(x_pix).astype(int), np.round(y_pix).astype(int)
            x_int, y_int = x_pix.astype(int), y_pix.astype(int)
            i, j = ax + x_int[mask], ay + y_int[mask]
            slice_array[i, j] = slice_photons[mask]

        else:
            # If astrometric precision is not that important and everything
            # has been oversampled, use this section.
            #  - ax, ay are the pixel coordinates of the image centre

            x_int, y_int = x_pix.astype(int), y_pix.astype(int)
            i, j = ax + x_int[mask], ay + y_int[mask]
            for ii, jj, ph in zip(i, j, slice_photons[mask]):
                slice_array[ii, jj] += ph

            try:
                slice_array = convolve_fft(slice_array, psf.array,
                                           allow_huge=True)
            except ValueError:
                slice_array = convolve(slice_array, psf.array)

        return slice_array


    def photons_in_range(self, lam_min=None, lam_max=None, min_bins=10,
                         mask=None):
        ## Argument 'mask' is unused (OC)
        """
        Calculate how many photons for each source exist in the wavelength range
        defined by lam_min and lam_max.

        Parameters
        ----------
        lam_min, lam_max : float, optional
            [um] integrate photons between these two limits. If both are `None`,
            limits are set at lam[0], lam[-1] for the source's wavelength range
        min_bins : float, optional
            the minimum number of spectral bins counted per layer
        """
        if lam_min is None:
            lam_min = self.lam[0]
        if lam_max is None:
            lam_max = self.lam[-1]

        # Check if the slice limits are within the spectrum wavelength range
        if lam_min > self.lam[-1] or lam_max < self.lam[0]:
            print((lam_min, lam_max), (self.lam[0], self.lam[-1]))
            raise ValueError("lam_min or lam_max outside wavelength range" + \
                                                                "of spectra")



        # find the closest indices i0, i1 which match the limits
        #x0, x1 = np.abs(self.lam - lam_min), np.abs(self.lam - lam_max)
        #i0 = np.where(x0 == np.min(x0))[0][0]
        #i1 = np.where(x1 == np.min(x1))[0][0]
        i0 = np.argmin(np.abs(self.lam - lam_min))
        i1 = np.argmin(np.abs(self.lam - lam_max))
        if self.lam[i0] > lam_min and i0 > 0:
            i0 -= 1
        if self.lam[i1] < lam_max and i1 < len(self.lam):
            i1 += 1

        # If there are less than min_bins between i0 and i1, then interpolate
        if i1 - i0 < min_bins:
            lam_zoom = np.linspace(lam_min, lam_max, min_bins)
            spec_zoom = np.zeros((self.spectra.shape[0], len(lam_zoom)))
            for i in range(len(self.spectra)):
                spec_zoom[i, :] = np.interp(lam_zoom, self.lam[i0:i1],
                                            self.spectra[i, i0:i1])
            spec_photons = np.sum(spec_zoom, axis=1)
        else:
            spec_photons = np.sum(self.spectra[:, i0:i1], axis=1)

        slice_photons = spec_photons[self.ref] * self.weight
        return slice_photons


    def apply_transmission_curve(self, transmission_curve):
        """
        Apply the values from a TransmissionCurve object to self.spectra

        Keywords:
        - transmission_curve: The TransmissionCurve to be applied
        """
        tc = deepcopy(transmission_curve)
        tc.resample(self.lam)
        self.spectra *= tc.val


    def _convert_to_photons(self):
        """
        convert the spectra to photons/(s m2)
        if [arcsec] are in the units, we want to find the photons per pixel
        if [um] are in the units, we want to find the photons per wavelength bin

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Come back and put in other energy units like Jy, mag, ergs !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        """
        self.units = u.Unit(self.params["units"])
        bases = self.units.bases

        #print(self.units, bases, self.area)

        factor = 1.
        if u.s not in bases:
            factor /= (self.exptime*u.s)
        if u.m not in bases:
            factor /= (self.area   *u.m**2)
        if u.micron in bases:
            factor *= (self.lam_res*u.um)
        if u.arcsec in bases:
            factor *= (self.pix_res*u.arcsec)**2

        #print((factor*self.units).unit, factor)

        self.units = (factor*self.units).unit
        self.spectra *= factor


    def _from_cube(self, filename):
        """
        Make a Source object from a cube in memory or a FITS cube on disk
        """
        if isinstance(filename, str) and os.path.exists(filename):
            hdr = fits.getheader(filename)
            cube = fits.getdata(filename)
        else:
            raise ValueError(filename+" doesn't exist")

        lam_res = hdr["CDELT3"]
        lam_min = hdr["CRVAL3"] - hdr["CRPIX3"] * lam_res
        lam_max = lam_min + hdr["NAXIS3"] * lam_res

        flux_map = np.sum(cube, axis=0).astype(dtype=np.float32)
        x, y = np.where(flux_map != 0)

        self.lam = np.linspace(lam_min, lam_max, hdr["NAXIS3"])
        ## TODO: ipt is undefined (OC)
        self.spectra = np.swapaxes(ipt[:, x, y], 0, 1)
        self.x = x
        self.y = y
        self.ref = np.arange(len(x))
        self.weight = np.ones(len(x))

        if "BUNIT" in hdr.keys():
            self.params["units"] = u.Unit(hdr["BUNIT"])
        if "EXPTIME" in hdr.keys():
            self.params["exptime"] = hdr["EXPTIME"]
        if "AREA"   in hdr.keys():
            self.params["area"] = hdr["AREA"]
        if "CDELT1" in hdr.keys():
            self.params["pix_res"] = hdr["CDELT1"]
        if "CUNIT1" in hdr.keys():
            self.params["pix_unit"] = hdr["CUNIT1"]
        self.lam_res = lam_res

        self._convert_to_photons()

    def _from_arrays(self, lam, spectra, x, y, ref, weight=None):
        """
        Make a Source object from a series of lists
        - x,y : [arcsec]

        """
        self.lam = lam
        self.spectra = spectra
        self.x = x
        self.y = y
        self.ref = ref
        self.weight = weight   if weight is not None   else np.array([1]*len(x))
        self.lam_res = np.median(lam[1:] - lam[:-1])

        if len(spectra.shape) == 1:
            self.spectra = np.array((spectra, spectra))

        self._convert_to_photons()


    def read(self, filename):
        """
        Read in a previously saved Source FITS file
        """
        ipt = fits.open(filename)
        dat0 = ipt[0].data
        hdr0 = ipt[0].header
        dat1 = ipt[1].data
        hdr1 = ipt[1].header
        ipt.close()

        self.x = dat0[0, :]
        self.y = dat0[1, :]
        self.ref = dat0[2, :]
        self.weight = dat0[3, :]

        lam_min, lam_max = hdr1["LAM_MIN"], hdr1["LAM_MAX"]
        self.lam_res = hdr1["LAM_RES"]
        self.lam = np.linspace(lam_min, lam_max, hdr1["NAXIS1"])
        self.spectra = dat1

        if "BUNIT" in hdr0.keys():
            self.params["units"] = u.Unit(hdr0["BUNIT"])
        if "EXPTIME" in hdr0.keys():
            self.params["exptime"] = hdr0["EXPTIME"]
        if "AREA"   in hdr0.keys():
            self.params["area"] = hdr0["AREA"]
        if "CDELT1" in hdr0.keys():
            self.params["pix_res"] = hdr0["CDELT1"]
        if "CUNIT1" in hdr0.keys():
            self.params["pix_unit"] = u.Unit(hdr0["CUNIT1"])
        self.lam_res = hdr1["LAM_RES"]

        self._convert_to_photons()

    def write(self, filename):
        """
        Write the current Source object out to a FITS file

        Parameters
        ----------
        filename : str
            where to save the FITS file

        Notes
        -----
        Just a place holder so that I know what's going on with the input table
        * The fist extension [0] contains an "image" of size 4 x N where N is the
        amount of sources. The 4 columns are x, y, ref, weight.
        * The second extension [1] contains an "image" with the spectra of each
        source. The image is M x len(spectrum), where M is the number of unique
        spectra in the source list. M = max(ref) - 1
        """

        # hdr = fits.getheader("../../../PreSim/Input_cubes/GC2.fits")
        # ipt = fits.getdata("../../../PreSim/Input_cubes/GC2.fits")
        # flux_map = np.sum(ipt, axis=0).astype(dtype=np.float32)
        # x,y = np.where(flux_map != 0)
        # ref = np.arange(len(x))
        # weight = np.ones(len(x))
        # spectra = np.swapaxes(ipt[:,x,y], 0, 1)
        # lam = np.linspace(0.2,2.5,231)

        xyHDU = fits.PrimaryHDU(np.array((self.x, self.y, self.ref, self.weight)))
        xyHDU.header["X_COL"] = "1"
        xyHDU.header["Y_COL"] = "2"
        xyHDU.header["REF_COL"] = "3"
        xyHDU.header["W_COL"] = "4"

        xyHDU.header["BUNIT"] = self.units.to_string()
        xyHDU.header["EXPTIME"] = self.params["exptime"]
        xyHDU.header["AREA"] = self.params["area"]
        xyHDU.header["CDELT1"] = self.params["pix_res"]
        xyHDU.header["CDELT2"] = self.params["pix_res"]
        xyHDU.header["CUNIT1"] = self.params["pix_unit"]
        xyHDU.header["CUNIT2"] = self.params["pix_unit"]

        xyHDU.header["SIMCADO"] = "SOURCE"

        specHDU = fits.ImageHDU(self.spectra)
        specHDU.header["CRVAL1"] = self.lam[0]
        specHDU.header["CRPIX1"] = 0
        specHDU.header["CDELT1"] = (self.lam_res, "[um] Spectral resolution")
        specHDU.header["LAM_MIN"] = (self.lam[0], "[um] Minimum wavelength")
        specHDU.header["LAM_MAX"] = (self.lam[-1], "[um] Maximum wavelength")
        specHDU.header["LAM_RES"] = (self.lam_res, "[um] Spectral resolution")

        hdu = fits.HDUList([xyHDU, specHDU])
        hdu.writeto(filename, clobber=True, checksum=True)


    def __str__(self):
        return "A photon source object"

    def __array__(self):
        if self.array is None:
            return np.zeros((self.naxis1, self.naxis2))
        else:
            return self.array

    def __getitem__(self, i):
        return (self.x[i], self.y[i],
                self.spectra[self.ref[i], :] * self.weight[i])

    def __mul__(self, x):
        newsrc = deepcopy(self)
        if isinstance(x, (sc.TransmissionCurve, sc.EmissionCurve,
                          sc.UnityCurve, sc.BlackbodyCurve)):
            newsrc.apply_transmission_curve(x)
        else:
            newsrc.array *= x
        return newsrc

    def __add__(self, x):
        newsrc = deepcopy(self)
        if isinstance(x, Source):
            if self.units != x.units:
                raise ValueError("units are not compatible: " + \
                                  str(self.units) + ", " + str(x.units))

            newsrc.lam = self.lam
            newsrc.spectra = list(self.spectra)
            for spec in x.spectra:
                tmp = np.interp(self.lam, x.lam, spec)
                newsrc.spectra += [tmp]
            newsrc.spectra = np.asarray(newsrc.spectra)

            newsrc.x = np.array((list(self.x) + list(x.x)))
            newsrc.y = np.array((list(self.y) + list(x.y)))
            newsrc.ref = np.array((list(self.ref) + list(x.ref + x.spectra.shape[0])))
            newsrc.weight = np.array((list(self.weight) + list(x.weight)))

            newsrc.x_orig = deepcopy(newsrc.x)
            newsrc.y_orig = deepcopy(newsrc.y)

        else:
            newsrc.array += x
        return newsrc

    def __sub__(self, x):
        newsrc = deepcopy(self)
        newsrc.array -= x
        return newsrc

    def __rmul__(self, x):
        return self.__mul__(x)

    def __radd__(self, x):
        return self.__add__(x)

    def __rsub__(self, x):
        return self.__mul__(-1) + x

    def __imul__(self, x):
        return self.__mul__(x)

    def __iadd__(self, x):
        return self.__add__(x)

    def __isub__(self, x):
        return self.__sub__(x)


##############################################################################




def _get_stellar_properties(spec_type, cat=None, verbose=False):
    """
    Returns an astropy.Table with the list of properties for the star in
    `spec_type`

    Parameters
    ----------
    spec_type : str, list
        The single or list of spectral types

    Returns
    -------
    props : astropy.Table or list of astropy.Tables with stellar paramters

    """

    if cat is None:
        cat = ioascii.read(os.path.join(__pkg_dir__, "data",
                                        "EC_all_stars.csv"))

    if isinstance(spec_type, (list, tuple)):
        return [_get_stellar_properties(i, cat) for i in spec_type]
    else:
        # Check if stellar type is in cat; if not look for the next
        # type in the sequence that is and assign its values
        spt, cls, lum = spec_type[0], int(spec_type[1]), spec_type[2:]
        for i in range(10):
            if cls > 9:
                cls = 0
                spt = "OBAFGKMLT"["OBAFGKMLT".index(spt)+1]

            startype = spt+str(cls)+lum # was 'star', redefined function star()
            cls += 1

            if startype in cat["Stellar_Type"]:
                break

        else:   # for loop did not find anything
            raise ValueError(spec_type+" doesn't exist in the database")

        n = np.where(cat["Stellar_Type"] == startype.upper())[0][0]
        if verbose:
            print("Returning properties for", startype)

        return cat[n]


def _get_stellar_mass(spec_type):
    """
    Returns a single (or list of) float(s) with the stellar mass(es) in units
    of Msol

    Parameters
    ----------
    spec_type : str, list
        The single or list of spectral types

    Returns
    -------
    mass : float, list

    """

    props = _get_stellar_properties(spec_type)

    if isinstance(props, (list, tuple)):
        return [prop["Mass"] for prop in props]
    else:
        return props["Mass"]


def _get_stellar_Mv(spec_type):
    """
    Returns a single (or list of) float(s) with the V-band absolute magnitude(s)

    Parameters
    ----------
    spec_type : str, list
        The single or list of spectral types

    Returns
    -------
    Mv : float, list

    """

    props = _get_stellar_properties(spec_type)

    if isinstance(props, (list, tuple)):
        return [prop["Mv"] for prop in props]
    else:
        return props["Mv"]


def _get_pickles_curve(spec_type, cat=None, verbose=False):
    """
    Returns the emission curve for a single or list of `spec_type`, normalised
    to 5556A

    Parameters
    ----------
    spec_type : str, list
        The single (or list) of spectral types (i.e. "A0V" or ["K5III", "B5I"])

    Returns
    -------
    lam : np.array
        a single np.ndarray for the wavelength bins of the spectrum,
    val : np.array (list)
        a (list of) np.ndarray for the emission curve of the spectral type(s)
        relative to the flux at 5556A

    References
    ----------
    Pickles 1998 - DOI: 10.1086/316197

    """
    if cat is None:
        cat = fits.getdata(os.path.join(__pkg_dir__, "data", "EC_pickles.fits"))

    if isinstance(spec_type, (list, tuple)):
        return cat["lam"], [_get_pickles_curve(i, cat)[1] for i in spec_type]
    else:
        # split the spectral type into 3 components and generalise for Pickles
        spt, cls, lum = spec_type[0], int(spec_type[1]), spec_type[2:]
        if lum.upper() == "I":
            lum = "Ia"
        elif lum.upper() == "II":
            lum = "III"
        elif "V" in lum.upper():
            lum = "V"


        for i in range(10):  # TODO: What does this loop do? (OC)
            if cls > 9:
                cls = 0
                spt = "OBAFGKMLT"["OBAFGKMLT".index(spt)+1]
            startype = spt+str(cls)+lum
            cls += 1

            if startype in cat.columns.names:
                break

        if spec_type != startype and verbose:
            print(spec_type, "isn't in Pickles. Returned", startype)

        try:
            lam, spec = cat["lam"], cat[startype]
        except KeyError:      # Correct? This shouldn't use error handling.
            lam, spec = cat["lam"], cat["M9III"]
        return lam, spec


def _scale_pickles_to_photons(spec_type, mag=0):
    """
    Pull in a spectrum from the Pickles library and scale the photon flux to a
    V=0 star

    Parameters
    ----------
    spec_type : str
        A spectral types (i.e. "A0V")

    Notes
    -----
    - Vega has a 5556 flux of between 950 and 1000 ph/s/cm2/A. The pickles
    resolution is 5 Ang.
    - Therefore the flux at 5555 should be 5 * 1000 * 10^(-0.4*Mv) ph/s/cm2/bin
    - Pickles catalogue is in units of Flambda [erg/s/cm2/A]
    - Ergo we need to divide the pickels values by lam/0.5556[nm], then rescale
    Regarding the number of photons in the 1 Ang bin at 5556 Ang
    - Bohlin (2014) says F(5556)=3.44×10−9 erg cm−2 s−1 A−1
    - Values range from 3.39 to 3.46 with the majority in range 3.44 to 3.46.
      Bohlin recommends 3.44
    - This results in a photon flux of 962 ph cm-2 s-1 A-1 at 5556 Ang
    """

    if isinstance(spec_type, (list, tuple, np.ndarray)):
        if isinstance(mag, (list, tuple, np.ndarray)):
            mag = list(mag)
        else:
            mag = [mag]*len(spec_type)
    else:
        mag = [mag]

    Mv = _get_stellar_Mv(spec_type)
    if not hasattr(Mv, "__len__"):
        Mv = [Mv]
    lam, ec = _get_pickles_curve(spec_type)
    dlam = (lam[1:] - lam[:-1])
    dlam = np.append(dlam, dlam[-1])

    lam *= 1E-4         # convert to um from Ang

    # Use Bohlin (2014) to determine the photon flux of a mag 0 A0V star
    # at 5556 Ang
    F = 3.44E-9 * u.erg / (u.cm**2 * u.s * u.AA)
    E = c.c*c.h/(5556*u.AA)
    ph0 = (F/E).to(1/(u.s * u.cm**2 * u.AA)).value

    # 5 Ang/bin * ~962 ph/s * (abs mag + apparent mag)
    # TODO: Wouldn't
    #   ph_factor = dlam * ph0 * 10**(-0.4 * (Mv + mag))
    # work?  (OC)
    ph_factor = []
    for i in range(len(mag)):
        tmp = dlam * ph0 * 10**(-0.4*(Mv[i] + mag[i]))
        ph_factor += [tmp]

    # take care of the conversion to ph/s/m2 by multiplying by 1E4
    # TODO: The original type(ec) == (list, tuple) is wrong (should be 'in')
    #   However, correcting it (using idiomatic isinstance) breaks the code!
    #   There must be a bug.
    # Correct code:
    # if isinstance(ec, (list, tuple)):
    #     for i in range(len(ec)):
    if type(ec) == (list, tuple):
        for i in len(range(ec)):
            ec[i] *= (lam/0.5556) * ph_factor[i] * 1E4
    else:
        ec *= (lam/0.5556) * ph_factor[0] * 1E4

    return lam, ec


def zero_magnitude_photon_flux(filter_name, area=1):
    """
    Return the number of photons for a m=0 star for a certain filter

    Parameters
    ----------
    filter_name : str
        the name of the broadband filter - UBVRIYzJHKKs
    area : float
        [m2] Collecting area of main mirror
    Notes
    -----
    units in [ph/s/m2]
    """

    if filter_name not in "UBVRIYzJHKKs":
        raise ValueError("Filter name must be one of UBVRIYzJHKKs: "+filter_name)

    lam, vega = _scale_pickles_to_photons("A0V", mag=-0.58)

    vraw = ioascii.read(os.path.join(__pkg_dir__, "data",
                                     "TC_filter_"+filter_name+".dat"))
    vlam = vraw[vraw.colnames[0]]
    vval = vraw[vraw.colnames[1]]
    filt = np.interp(lam, vlam, vval)

    n_ph = np.sum(vega*filt) * area

    #print("units in [ph/s/m2]")
    return n_ph


def value_at_lambda(lam_i, lam, val, return_index=False):
    """
    Return the value at a certain wavelength - i.e. val[lam] = x

    Parameters
    ----------
    lam_i : float
        the wavelength of interest
    lam : np.ndarray
        an array of wavelengths
    val : np.ndarray
        an array of values
    return_index : bool, optional
        If True, the index of the wavelength of interest is returned
        Default is False
    """

    i0 = np.where(lam <= lam_i)[0][-1]
    i1 = np.where(lam > lam_i)[0][0]

    lam_x = np.array([lam[i0], lam_i, lam[i1]])
    val_i = np.interp(lam_x, lam, val)

    if return_index:
        return i0
    else:
        return val_i[1]


def SED(spec_type, filter_name="V", magnitude=0.):
    """
    Return a scaled SED for a X type star at magnitude Y star in band Z

    Parameters
    ----------
    spec_type : str
        The spectral type of the star
    filter_name : str, optional
        [UBVRIYzJHKKs] Any braodband filter in the Vis+NIR regime.
        Default is "V"
    magnitude : float, optional
        Apparent magnitude of the star. Default is 0.

    Returns
    -------
    lam : np.ndarray
        [um] The centre of each 5 Ang bin along the spectral axis
    val : np.ndarray
        [ph/s/m2/bin] The photon flux of the star in each bin

    Notes
    -----
    Original flux units are in [ph/s/m2/AA], so we multiply the flux by 5 to
    get [ph/s/m2/bin]. Therefore divide by 5*1E4 if you need the flux in
    [ph/s/cm2/Angstrom]

    """

    if filter_name not in "UBVRIYzJHKKs":
        raise ValueError("Filter name must be one of UBVRIYzJHKKs: "+filter_name)

    if isinstance(magnitude, (list, tuple)):
        magnitude = np.asarray(magnitude)

    flux_0 = zero_magnitude_photon_flux(filter_name)
    flux = flux_0 * 10**(-0.4 * magnitude)

    lam, starflux = _scale_pickles_to_photons(spec_type)
    # was 'star' but that redefines function star()

    vraw = ioascii.read(os.path.join(__pkg_dir__, "data",
                                     "TC_filter_"+filter_name+".dat"))
    vlam = vraw[vraw.colnames[0]]
    vval = vraw[vraw.colnames[1]]
    filt = np.interp(lam, vlam, vval)

    n_ph = np.sum(starflux  * filt)

    scale_factor = flux / n_ph
    #print("scale_factor, flux, n_ph, flux_0 [ph/s/m2/bin]")
    #print(scale_factor, flux, n_ph, flux_0, filter_name)

    return lam, (scale_factor * starflux.transpose()).transpose()


def empty_sky():
    """
    Returns an empty source so that instrumental fluxes can be simulated
    """
    return Source(lam=np.linspace(0.3, 2.5, 221),
                  spectra=np.zeros((2, 221)),
                  x=[0], y=[0], ref=[0], weight=[0])


def star_grid(n, mag_min, mag_max, filter_name="K", separation=1, area=1,
              spec_type="A0V"):
    """
    Creates a square grid of A0V stars at equal magnitude intervals

    Parameters
    ----------
    n : float
        the number of stars in the grid
    mag_min, mag_max : float
        [vega mag] the minimum (brightest) and maximum (faintest) magnitudes for
        stars in the grid
    filter_name : str
        broadband filter - B V R I Y z J H K Ks
    separation : float, optional
        [arcsec] an average speration between the stars in the grid can be
        specified. Default is 1 arcsec
    area : float, optional
        [m2] collecting area of primary mirror
    spec_type : str, optional
        the spectral type of the star, e.g. "A0V", "G5III"

    Returns
    -------
    source : `simcado.Source`

    Notes
    -----
    The units of the A0V spectrum in `source` are [ph/s/bin] or [ph/s/m2/bin]
    depending on if area is given.
    The weight values are the scaling factors to bring a V=0 A0V spectrum down
    to the required magnitude for each star.

    """

    if isinstance(mag_min, (list, tuple, np.ndarray)):
        mags = np.asarray(mag_min)
    else:
        if mag_min < mag_max:
            mags = np.linspace(mag_min, mag_max, n)
        elif mag_min < mag_max:
            mags = np.linspace(mag_max, mag_min, n)
        elif mag_min == mag_max:
            mags = np.ones(n) * mag_min

    side_len = int(np.sqrt(n)) + (np.sqrt(n) % 1 > 0)

    x = separation * (np.arange(n) % side_len - (side_len - 1) / 2)
    y = separation * (np.arange(n)// side_len - (side_len - 1) / 2)

    lam, spec = SED(spec_type, filter_name=filter_name, magnitude=0)
    if isinstance(spec_type, (list, tuple)):
        ref = np.arange(len(spec_type))
    else:
        ref = np.zeros((n))
    weight = 10**(-0.4*mags) * area

    if area == 1:
        units = "ph/s/m2"
    else:
        units = "ph/s"

    src = Source(lam=lam, spectra=spec,
                 x=x, y=y,
                 ref=ref, weight=weight,
                 units=units,
                 area=area)

    return src


def star(mag, filter_name="K", spec_type="A0V", position=(0, 0)):
    """
    Creates a simcado.Source object for a star with a given magnitude

    Parameters
    ----------
    mag : float
        magnitude of star
    filter_name : str
        filter in which the magnitude is given
    spec_type : str, optional
        the spectral type of the star, e.g. "A0V", "G5III"

    Returns
    -------
    source : `simcado.Source`
    """
    thestar = star_grid(1, mag, mag, filter_name, spec_type=spec_type)
    thestar.x, thestar.y = [position[0]], [position[1]]
    return thestar


def stars(mags, x, y, filter_name="K", spec_types="A0V"):
    """
    Creates a simcado.Source object for a bunch of stars.

    Parameters
    ----------
    mags : array
         magnitudes of the stars
    x, y : arrays
         x and y coordinates of the stars
    filter_name : str
         filter in which the magnitudes are given
    spec_type : str or array of strings (optional)
        the spectral type(s) of the stars, e.g. "A0V", "G5III"

    Returns
    -------
    source : `simcado.Source`
    """
    if isinstance(spec_types, (tuple, list)) and len(mags) != len(spec_types):
        raise ValueError("len(mags) != len(spec_types)")

    thestars = star_grid(len(mags), mags, mags, filter_name,
                         spec_type=spec_types)
    thestars.x, thestars.y = x, y
    return thestars


def source_1E4_Msun_cluster(distance=50000, half_light_radius=1):
    """
    Generate a source object for a 10^4 solar mass cluster

    Parameters
    ----------
    distance : float
        [pc] distance to the cluster
    half_light_radius : float
        [pc] half light radius of the cluster

    Returns
    -------
    src : simcado.Source

    """
    # IMF is a realisation of stellar masses drawn from an initial mass
    # function (TODO: which one?) summing to 1e4 M_sol.
    fname = os.path.join(__pkg_dir__, "data", "IMF_1E4.dat")
    imf = np.loadtxt(fname)

    # Assign stellar types to the masses in imf using list of average
    # main-sequence star masses:
    stel_type = [i + str(j) + "V" for i in "OBAFGKM" for j in range(10)]
    mass = _get_stellar_mass(stel_type)
    ref = utils.nearest(mass, imf)
    thestars = [stel_type[i] for i in ref] # was stars, redefined function name

    # assign absolute magnitudes to stellar types in cluster
    unique_ref = np.unique(ref)
    unique_type = [stel_type[i] for i in unique_ref]
    unique_Mv = _get_stellar_Mv(unique_type)

    # Mv_dict = {i : float(str(j)[:6]) for i, j in zip(unique_type, unique_Mv)}
    ref_dict = {i : j for i, j in zip(unique_type, np.arange(len(unique_type)))}

    # find spectra for the stellar types in cluster
    lam, spectra = _scale_pickles_to_photons(unique_type)

    # this one connects the stars to one of the unique spectra
    stars_spec_ref = [ref_dict[i] for i in thestars]

    # absolute mag + distance modulus
    m = np.array([unique_Mv[i] for i in stars_spec_ref])
    m += 5 * np.log10(distance) - 5

    # set the weighting
    weight = 10**(-0.4*m)

    # draw positions of stars: cluster has Gaussian profile
    distance *= u.pc
    half_light_radius *= u.pc
    hwhm = (half_light_radius/distance*u.rad).to(u.arcsec).value
    sig = hwhm / np.sqrt(2 * np.log(2))

    x = np.random.normal(0, sig, len(imf))
    y = np.random.normal(0, sig, len(imf))

    src = Source(lam=lam, spectra=spectra, x=x, y=y, ref=stars_spec_ref,
                 weight=weight, units="ph/s/m2")

    return src



def source_from_image(images, lam, spectra, pix_res, oversample=1,
                      units="ph/s/m2", flux_threshold=0,
                      center_pixel_offset=(0, 0)):
    """
    Create a Source object from an image or a list of images.

    Parameters
    ----------
    images : np.ndarray, list
        A single or list of np.ndarrays describing where the flux is coming from.
        The spectrum for each pixel in the image is weighted by the pixel value.
    lam : np.ndarray
        An array contains the centres of the wavelength bins for the spectra
    spectra : np.ndarray
        A (n,m) array with n spectra, each with m bins
    pix_res : float
        [arcsec] The pixel scale of the images in arcseconds
    oversample : int
        The factor with which to oversample the image. Each image pixel is split
        into (oversample)^2 individual point sources.
    units : str, optional
        The energy units of the spectra. Default is [ph/s/m2]
    flux_threshold : float, optional
        If there is noise in the image, set threshold to the noise limit so that
        only real photon sources are extracted. Default is 0.
    center_pixel_offset : (int, int)
        [pixel] If the central pixel is offset from the centre of the image, add
        this offset to (x,y) coordinates.

    Returns
    -------
    src : source.Source


    """

    if isinstance(images, (list, tuple)):
        srcs = [source_from_image(images[i], lam, spectra[i, :], pix_res,
                                  oversample, units, flux_threshold,
                                  center_pixel_offset)
                for i in range(len(images))]
        src = srcs[0]
        for i in range(1, len(images)):
            src += srcs[i]

    else:
        if not isinstance(oversample, int):
            raise ValueError("Oversample must be of type 'int'")

        if isinstance(images, str) and images.split(".")[-1].lower() == "fits":
            images = fits.getdata(images)

        im = images
        x_cen, y_cen = np.array(im.shape) / 2 + np.array(center_pixel_offset)
        x_i, y_i = np.where(im > flux_threshold)

        x = (x_i - x_cen) * pix_res
        y = (y_i - y_cen) * pix_res
        weight = im[x_i, y_i]

        i = oversample
        oset = np.linspace(-0.5, 0.5, 2*i+1)[1:2*i:2] * pix_res

        x_list, y_list, w_list = [], [], []
        for i in oset:
            for j in oset:
                x_list += (x + i).tolist()
                y_list += (y + j).tolist()
                w_list += (weight / oversample**2).tolist()
        x, y, weight = np.array(x_list), np.array(y_list), np.array(w_list)

        ref = np.zeros(len(x))

        src = Source(lam=lam, spectra=spectra, x=x, y=y, ref=ref, weight=weight,
                     units=units)

        return src


def scale_spectrum(lam, spec, mag, filter_name="K", return_ec=False):
    """
    Scale a spectrum to be a certain magnitude

    Parameters
    ----------
    lam : np.ndarray
        [um] The wavelength bins for spectrum
    spec : np.ndarray
        The spectrum to be scaled into [ph/s/m2] for the given broadband filter
    mag : float
        magnitude of the source
    filter_name : str, optional
        broadband filter in the Vis/NIR range UBVRIzYJHKKs. Default is "K"
    return_ec : bool, optional
        If True, a simcado.spectral.EmissionCurve object is returned.
        Default is False

    Returns
    -------
    lam : np.ndarray
        [um] The centres of the wavelength bins for the new spectrum
    spec : np.array
        [ph/s/m2] The spectrum scaled to the specified magnitude

    """

    ideal_ph = zero_magnitude_photon_flux(filter_name) * 10**(-0.4 * mag)

    from simcado.spectral import EmissionCurve
    from simcado.optics import filter_curve

    curve = EmissionCurve(lam=lam, val=spec, area=1, units="ph/s/m2")
    curve.resample(0.001)

    filt = filter_curve(filter_name)

    tmp = curve * filt
    obs_ph = tmp.photons_in_range()

    scale_factor = ideal_ph / obs_ph
    curve *= scale_factor

    if return_ec:
        return curve
    else:
        return curve.lam, curve.val


def scale_spectrum_sb(lam, spec, mag_per_arcsec, pix_res=0.004, filter_name="K",
                      return_ec=False):
    """
    Scale a spectrum to be a certain magnitude per arcsec

    Parameters
    ----------
    lam : np.ndarray
        [um] The wavelength bins for spectrum
    spec : np.ndarray
        The spectrum to be scaled into [ph/s/m2] for the given broadband filter
    mag_per_arcsec : float
        [mag] surface brightness of the source
    pix_res : float
        [arcsec] the pixel resolution
    filter_name : str, optional
        broadband filter in the Vis/NIR range UBVRIzYJHKKs. Default is "K"
    return_ec : bool, optional
        If True, a simcado.spectral.EmissionCurve object is returned.
        Default is False

    Returns
    -------
    lam : np.ndarray
        [um] The centres of the wavelength bins for the new spectrum
    spec : np.array
        [ph/s/m2/arcsec] The spectrum scaled to the specified magnitude

    """

    if return_ec:
        curve = scale_spectrum(lam, spec, mag_per_arcsec, filter_name,
                               return_ec)
        curve.val *= pix_res**2
        curve.params["pix_res"] = pix_res
        return curve
    else:
        lam, spec = scale_spectrum(lam, spec, mag_per_arcsec, filter_name,
                                   return_ec)
        spec *= pix_res**2
        return lam, spec


def flat_spectrum(mag, filter_name="K", return_ec=False):
    """
    Return a flat spectrum scaled to a certain magnitude

    Parameters
    ----------
    mag : float
        [mag] magnitude of the source
    filter_name : str, optional
        broadband filter in the Vis/NIR range UBVRIzYJHKKs. Default is "K"
    return_ec : bool, optional
        If True, a simcado.spectral.EmissionCurve object is returned.
        Default is False

    Returns
    -------
    lam : np.ndarray
        [um] The centres of the wavelength bins for the new spectrum
    spec : np.array
        [ph/s/m2/arcsec] The spectrum scaled to the specified magnitude

    """
    lam = np.arange(0.3, 2.51, 0.01)
    spec = np.ones(len(lam))

    if return_ec:     # TODO: mag_per_arcsec undefined? (OC)
        curve = scale_spectrum(lam, spec, mag_per_arcsec, filter_name,
                               return_ec)
        return curve
    else:
        lam, spec = scale_spectrum(lam, spec, mag_per_arcsec, filter_name,
                                   return_ec)
        return lam, spec


def flat_spectrum_sb(mag_per_arcsec, filter_name="K", pix_res=0.004,
                     return_ec=False):
    """
    Return a flat spectrum for a certain magnitude per arcsec

    Parameters
    ----------
    mag : float
        [mag] magnitude of the source
    filter_name : str, optional
        broadband filter in the Vis/NIR range UBVRIzYJHKKs. Default is "K"
    pix_res : float
        [arcsec] the pixel resolution. Default is 4mas (i.e. 0.004)
    return_ec : bool, optional
        If True, a simcado.spectral.EmissionCurve object is returned. Default is False

    Returns
    -------
    lam : np.ndarray
        [um] The centres of the wavelength bins for the new spectrum
    spec : np.array
        [ph/s/m2/arcsec] The spectrum scaled to the specified magnitude

    """
    lam = np.arange(0.3, 2.51, 0.01)
    spec = np.ones(len(lam))

    if return_ec:
        curve = scale_spectrum(lam, spec, mag_per_arcsec, filter_name, return_ec)
        curve.val *= pix_res**2
        return curve
    else:
        lam, spec = scale_spectrum(lam, spec, mag_per_arcsec, filter_name,
                                   return_ec)
        spec *= pix_res**2
        return lam, spec
