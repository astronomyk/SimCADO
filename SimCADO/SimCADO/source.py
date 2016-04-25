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
from copy import deepcopy
import warnings

import numpy as np
import scipy.ndimage.interpolation as spi

from astropy.io import fits
from astropy.convolution import convolve, convolve_fft
import astropy.units as u

try:
    import SimCADO.spatial as pe
except:
    import spatial as pe


__all__ = ["source", "Source"]


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
    - pix_res
    - exptime
    - area
    """

    def __init__(self, filename=None,
                lam=None, spectra=None, x=None, y=None, ref=None, weight=None,
                **kwargs):

        self.params = {"units"   :"ph/s",
                       "pix_unit":"arcsec",
                       "pix_res" :0.004,
                       "exptime" :1,
                       "area"    :1}
        self.params.update(kwargs)

        self.info = dict([])
        self.info['created'] = 'yes'
        self.info['description'] = "List of spectra and their positions"

        self.units = u.Unit(self.params["units"])
        self.pix_res = self.params["pix_res"]
        self.exptime = self.params["exptime"]
        self.area = self.params["area"]

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

        self.ref = self.ref.astype(int)
        self.x_orig = deepcopy(self.x)
        self.y_orig = deepcopy(self.y)


    def apply_optical_train(self, opt_train, chips, **kwargs):
        """

        Output array is in units of [ph/s/pixel]


        Parameters
        ==========
        opt_train : OpticalTrain object
        chips : int
            Default is 0

        """
        params = {"verbose"     :False,
                  "ATMO_BG_ON"  :"yes",
                  "INST_DEROT_PERFORMANCE"  :100,
                  "SCOPE_JITTER_FWHM"       :0,
                  "SCOPE_DRIFT_DISTANCE"    :0     }
        params.update(kwargs)

        # 0. Create a canvas onto which we splat the PSFed sources 
        # 1. Apply the master transmission curve to all the spectra
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
        # 4. Add the average number of atmo-bg and mirror-bb photons
        # 5. Apply the instrumental distortion

        if not hasattr(chips, "__len__"):
            chips = [chips]
        
        for chip_i in chips:
    
            # 0.
            image = None

            # 1.
            self.apply_transmission_curve(opt_train.tc_source)

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
                psf = opt_train.psf_source[i]

                if image == None:
                    image = self.image_in_range(psf, lam_min, lam_max,
                                                opt_train.chips[chip_i])
                else:
                    image += self.image_in_range(psf, lam_min, lam_max,
                                                 opt_train.chips[chip_i])

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

            # 4.
            image += (opt_train.n_ph_atmo + opt_train.n_ph_mirror)

            self.project_onto_chip(image, opt_train.chips[chip_i])

        ######################################
        # CAUTION WITH THE PSF NORMALISATION #
        ######################################


    def project_onto_chip(self, image, chip):

        scale_factor = self.pix_res / chip.pix_res
        chip_arr = spi.zoom(image, scale_factor, order=1)
        chip_arr *= np.sum(image) / np.sum(chip_arr)

        chip.add_signal(chip_arr)


    def image_in_range(self, psf, lam_min, lam_max, chip=None, **kwargs):
        """
        Find the sources that fall in the chip area and generate an image for
        the wavelength range [lam_min, lam_max)

        Output is in [ph/s/pixel]

        Parameters
        ==========
        psf : psf.PSF object
            The PSF that the sources will be convolved with
        lam_min, lam_max : float
            [um] the wavelength range relevant for the psf
        chip : detector.Chip object
            the chip that will be seeing this image. If chip is None, an image
            covering the entire Source region will be generated at a resolution
            of `pix_res`. This could take a long time!

        **kwargs
        ========
        sub_pixel : bool
            if sub-pixel accuracy is needed, each source is shifted individually.
            Default is False
        pix_res : float
            [arcsec] the field of view of each pixel. Default is 0.004 arcsec
        oversample : int
            the psf images will be oversampled to better conserve flux.
            Default is 1 (i.e. not oversampled)
        """

        params = {"pix_res"     :0.004,
                  "sub_pixel"   :False,
                  "oversample"  :1      }

        params.update(kwargs)

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
            x_cen, y_cen = (x_max + xmin) / 2, (y_max + ymin) / 2

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

        x_int, y_int = x_pix.astype(int), y_pix.astype(int)
        dx, dy = self.x - x_int, self.y - y_int

        if params["sub_pixel"]:
            # for each point source in the list, add a psf to the slice_array
            for i in range(len(slice_photons)):
                if not mask[i]:
                    continue

                psf_tmp = spi.shift(psf.array, (dx[i],dy[i]), order=1)
                x_pint, y_pint = x_int[i], y_int[i]

                # Find the slice borders for the array where the psf will go
                ax0 = np.max(np.array((x_pint - bx, [0]*len(x_pint))), axis=0)
                ax1 = np.min(np.array((x_pint + bx + 1,
                                 [slice_array.shape[0]]*len(x_pint))), axis=0)
                ay0 = np.max(np.array((y_pint - by, [0]*len(y_pint))), axis=0)
                ay1 = np.min(np.array((y_pint + by + 1,
                                 [slice_array.shape[1]]*len(y_pint))), axis=0)

                # the slice limits of the psf array are found by taking the
                # pixel distance from the x,y position to the slice limits
                # of the slice_array. This distance is subtracted from the
                # centre of the psf array.
                bx0 = bx - (x_pint - ax0)
                bx1 = bx + (ax1 - x_pint)
                by0 = by - (y_pint - ay0)
                by1 = by + (ay1 - y_pint)

                slice_array[ax0:ax1, ay0:ay1] = psf.array[bx0:bx1, by0:by1] \
                                        * slice_photons[i] * self.weights[i]

        else:
            # If astrometric precision is not that important and everything
            # has been oversampled, use this section.
            #  - ax, ay are the pixel coordinates of the image centre
            print(lam_min,lam_max)

            i, j = ax + x_int[mask], ay + y_int[mask]
            slice_array[i,j] = slice_photons[mask] * self.weight[mask]
            try:
                slice_array = convolve_fft(slice_array, psf.array, allow_huge=True)
            except:
                slice_array = convolve(slice_array, psf.array)

        return slice_array


    def photons_in_range(self, lam_min, lam_max, min_bins=10, mask=None):
        """
        Calculate how many photons for each source exist in the wavelength range
        defined by lam_min and lam_max.

        Keywords
        ========
        - lam_min, lam_max : float

        Optional keywords:
        - min_bins : float
            the minimum number of spectral bins counted per layer
        """
        # Check if the slice limits are within the spectrum wavelength range
        if lam_min > self.lam[-1] or lam_max < self.lam[0]:
            print((lam_min, lam_max), (self.lam[0], self.lam[-1]))
            raise ValueError("lam_min or lam_max outside wavelength range" + \
                                                                "of spectra")

        # find the closest indices i0, i1 which match the limits
        x0, x1 = np.abs(self.lam - lam_min), np.abs(self.lam - lam_max)
        i0 = np.where(x0 == np.min(x0))[0][0]
        i1 = np.where(x1 == np.min(x1))[0][0]
        if self.lam[i0] > lam_min and i0 > 0:
            i0 -= 1
        if self.lam[i1] < lam_max and i1 < len(self.lam):
            i1 += 1

        # If there are less than min_bins between i0 and i1, then interpolate
        if i1 - i0 < min_bins:
            lam_zoom  = np.linspace(lam_min, lam_max, min_bins)
            spec_zoom = np.zeros((self.spectra.shape[0], len(lam_zoom)))
            for i in range(len(self.spectra)):
                spec_zoom[i,:] = np.interp(lam_zoom, self.lam[i0:i1],
                                                        self.spectra[i,i0:i1])
            spec_photons = np.sum(spec_zoom, axis=1)
        else:
            spec_photons = np.sum(self.spectra[:,i0:i1], axis=1)

        slice_photons = spec_photons[self.ref]
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
        bases  = self.units.bases

        factor = 1.
        if u.s      not in bases: factor /= (self.exptime*u.s)
        if u.m      not in bases: factor /= (self.area   *u.m**2)
        if u.micron     in bases: factor *= (self.lam_res*u.um)
        if u.arcsec     in bases: factor *= (self.pix_res*u.arcsec)**2
        #print((factor*self.units).unit)

        self.units = (factor*self.units).unit
        self.spectra *= factor


    def _from_cube(self, filename, **kwargs):
        """
        Make a Source object from a cube in memory or a FITS cube on disk
        """
        if type(filename) == str and os.path.exists(filename):
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
        self.spectra = np.swapaxes(ipt[:,x,y], 0, 1)
        self.x = x
        self.y = y
        self.ref = np.arange(len(x))
        self.weight = np.ones(len(x))

        if "BUNIT"  in hdr.keys():  self.params["units"]   = u.Unit(hdr["BUNIT"])
        if "EXPTIME" in hdr.keys(): self.params["exptime"] = hdr["EXPTIME"]
        if "AREA"   in hdr.keys():  self.params["area"]    = hdr["AREA"]
        if "CDELT1" in hdr.keys():  self.params["pix_res"] = hdr["CDELT1"]
        if "CUNIT1" in hdr.keys():  self.params["pix_unit"] = hdr["CUNIT1"]
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

        self.x = dat0[0,:]
        self.y = dat0[1,:]
        self.ref = dat0[2,:]
        self.weight = dat0[3,:]

        lam_min, lam_max = hdr1["LAM_MIN"], hdr1["LAM_MAX"]
        self.lam_res     = hdr1["LAM_RES"]
        self.lam = np.linspace(lam_min, lam_max, hdr1["NAXIS1"])
        self.spectra = dat1

        if "BUNIT"  in hdr0.keys():     self.params["units"]   = u.Unit(hdr0["BUNIT"])
        if "EXPTIME" in hdr0.keys():    self.params["exptime"] = hdr0["EXPTIME"]
        if "AREA"   in hdr0.keys():     self.params["area"]    = hdr0["AREA"]
        if "CDELT1" in hdr0.keys():     self.params["pix_res"] = hdr0["CDELT1"]
        if "CUNIT1" in hdr0.keys():     self.params["pix_unit"] = u.Unit(hdr0["CUNIT1"])
        self.lam_res = hdr1["LAM_RES"]

        self._convert_to_photons()

    def write(self, filename):
        """
        Write the current Source object out to a FITS file

        Parameters:
        ===========
        filename:


        Just a place holder so that I know what's going on with the input table
        * The fist extension [0] contains an "image" of size 4 x N where N is the
        amount of sources. The 4 columns are x, y, ref, weight.
        * The second extension [1] contains an "image" with the spectra of each
        source. The image is M x len(spectrum), where M is the number of unique
        spectra in the source list. max(ref) = M - 1
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
        specHDU.header["CRVAL1"]  = self.lam[0]
        specHDU.header["CRPIX1"]  = 0
        specHDU.header["CDELT1"]  = (self.lam_res, "[um] Spectral resolution")
        specHDU.header["LAM_MIN"] = (self.lam[0], "[um] Minimum wavelength")
        specHDU.header["LAM_MAX"] = (self.lam[-1], "[um] Maximum wavelength")
        specHDU.header["LAM_RES"] = (self.lam_res, "[um] Spectral resolution")

        hdu = fits.HDUList([xyHDU, specHDU])
        hdu.writeto(filename, clobber=True)


    def __str__(self):
        return "A photon source object"

    def __array__(self):
        if self.array is None:
            return np.zeros((self.naxis1, self.naxis2))
        else:
            return self.array

    def __getitem__(self, i):
        return self.x[i], self.y[i], self.spectra[self.ref[i],:] * self.weight[i]

    def __mul__(self, x):
        newsrc = deepcopy(self)
        newsrc.array *= x
        return newsrc

    def __add__(self, x):
        newsrc = deepcopy(self)
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



class bloedsinn():
    pass
    
    
    

# class source(object):

    # def __init__(self, source, cmds):
        # """
        # Keywords:
        # - lam: [um] an array of size L with wavelengths for the spectra
        # - spectra: [photons] a 2D array of size (S,L) with a spectrum for S
                    # unique sources
        # - x, y: [pix] arrays of size N holding the pixel coordinate information
                # for N sources

        # Optional keywords:
        # - ref: [int] an array holding N references joining each of the N
                    # source pixel coordinates to one of the S unique spectra
                    # max(ref) < spectra.shape[0]
        # - weights: [float] an array of size N with weights for each source
        # """
        # self.cmds = cmds
        # # self.size = 512

        # self.info = dict([])
        # self.info['created'] = 'yes'
        # self.info['description'] = "List of spectra and their positions"

        # self.lam        = source.lam
        # self.spectra    = source.spectra
        # self.x_orig     = source.x
        # self.y_orig     = source.y
        # self.ref        = np.array(source.ref, dtype=int)
        # self.weight     = source.weight
        # self.x, self.y  = deepcopy(self.x_orig), deepcopy(self.y_orig)
        # self.src_params = source.params

        # # add a second dimension to self.spectra so that all the 2D calls work
        # if len(self.spectra.shape) == 1:
            # self.spectra = np.asarray([self.spectra]*2)

        # self.array = np.zeros((self.size, self.size), dtype=np.float32)

        # self.lam_bin_edges = cmds.lam_bin_edges
        # self.lam_bin_centers = cmds.lam_bin_centers


    # def _to_pixel_coords(self, ra_offset, dec_offset):
        # """
        # Return pixel coordinates
        # """
        # factor = 1.
        # if u.Unit(self.src_params["pix_unit"]) != u.arcsec:
            # factor = u.Unit(self.src_params["pix_unit"]).to(u.arcsec)

        # x_pix = factor * ra_offset / self.src_params["pix_res"]
        # y_pix = factor * dec_offset / self.src_params["pix_res"]

        # return x_pix, y_pix


    # def __str__(self):
        # return self.info['description']





    # def poissonify(self, arr=None):
        # """
        # Add a realisation of the poisson process to the array 'arr'.
        # If arr=None, the poisson realisation is applied to source.array

        # Optional keyword:
        # - arr:
        # """
        # if arr is None:
            # self.array = np.random.poisson(self.array).astype(np.float32)
        # else:
            # return np.random.poisson(arr).astype(np.float32)
