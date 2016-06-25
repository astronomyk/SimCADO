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

from astropy.io import fits, ascii
from astropy.convolution import convolve, convolve_fft
import astropy.units as u

try:
    import simcado.spatial as pe
    import simcado.spectral as sc
    __file__ = pe.__file__
except:
    import spatial as pe
    import spectral as sc
    __file__ = "./spatial.py"
    
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
        params = {"verbose"     :False,
                  "ATMO_BG_ON"  :"yes",
                  "INST_DEROT_PERFORMANCE"  :100,
                  "SCOPE_JITTER_FWHM"       :0,
                  "SCOPE_DRIFT_DISTANCE"    :0     }
        params.update(kwargs)

        
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

        if chips == None or str(chips).lower() == "all":
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
                psf = opt_train.psf_source[i]

                oversample = opt_train.cmds["SIM_OVERSAMPLING"]
                if image == None:
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

            # 5. 
            
            
            self.project_onto_chip(image, detector.chips[chip_i])

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
                  "oversample"  :1}

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
            slice_array[i,j] = slice_photons[mask]# * self.weight[mask]
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
        return self.x[i], self.y[i], self.spectra[self.ref[i],:] *self.weight[i]

    def __mul__(self, x):
        newsrc = deepcopy(self)
        
        if type(x) == sc.TransmissionCurve:
            newsrc.apply_transmission_curve(x)
        else:
            newsrc.array *= x
        return newsrc

    def __add__(self, x):
        newsrc = deepcopy(self)
        if type(x) == Source:
            if self.units != x.units:
                raise ValueError("units are not compatible: " + \
                                  str(self.units) + ", " + str(x.units))
                                  
            newsrc.lam = self.lam
            newsrc.spectra = self.spectra.tolist()
            for spec in x.spectra:
                tmp = np.interp(self.lam, x.lam, spec)
                newsrc.spectra += [tmp]
            newsrc.spectra = np.asarray(newsrc.spectra)
            
            newsrc.x = np.append(self.x, x.x)
            newsrc.y = np.append(self.y, x.y)
            newsrc.ref = np.append(self.ref, x.ref + x.spectra.shape[0])
            newsrc.weight = np.append(self.weight, x.weight)
            
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



import astropy.units as u
import astropy.constants as c


def get_stellar_properties(spec_type, cat=None, verbose=False):
    """
    Returns an astropy.Table with the list of properties for the star in 
    `spec_type`
    
    Parameters
    ----------
    spec_type : str, list
        The single or list of spectral types
        
    Returns
    -------
    props : astropy.Table or list of astropy.Tables with the paramters to the 
    
    """
    
    if cat is None:
        cat = ascii.read(os.path.join(__pkg_dir__, "data", "EC_all_stars.csv"))

    if type(spec_type) in (list, tuple):
        return [get_stellar_properties(i, cat) for i in spec_type]
    else:
        spt, cls, lum = spec_type[0], int(spec_type[1]), spec_type[2:]
        for i in range(10):
            if cls > 9:
                cls = 0
                spt = "OBAFGKMLT"["OBAFGKMLT".index(spt)+1]
            star = spt+str(cls)+lum
            cls += 1
            
            if star in cat["Stellar_Type"]: break 
        
        try:
            n = np.where(cat["Stellar_Type"] == star.upper())[0][0]
            if verbose: print("Returning properties for", star)
        except:
            raise ValueError(spec_type+" doesn't exist in the database")
        return cat[n]


def get_stellar_mass(spec_type):
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
    
    props = get_stellar_properties(spec_type)
    
    if type(props) in (list, tuple):
        return [prop["Mass"] for prop in props] 
    else:
        return props["Mass"]
    
    
def get_stellar_Mv(spec_type):
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
    
    props = get_stellar_properties(spec_type)
    
    if type(props) in (list, tuple):
        return [prop["Mv"] for prop in props] 
    else:
        return props["Mv"]
    
    
def get_pickles_curve(spec_type, cat=None, verbose=False):
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
        
    if type(spec_type) in (list, tuple):
        return cat["lam"], [get_pickles_curve(i, cat)[1] for i in spec_type]
    else:
        # split the spectral type into 3 components and generalise for Pickles
        spt, cls, lum = spec_type[0], int(spec_type[1]), spec_type[2:]
        if lum.upper() == "I": 
            lum = "Ia"
        elif lum.upper() == "II": 
            lum = "III"
        elif "V" in lum.upper(): 
            lum = "V"
        
        
        for i in range(10):
            if cls > 9:
                cls = 0
                spt = "OBAFGKMLT"["OBAFGKMLT".index(spt)+1]
            star = spt+str(cls)+lum
            cls += 1
            
            if star in cat.columns.names: break
                
        if spec_type != star and verbose:
            print(spec_type, "isn't in Pickles. Returned", star)
            
        return cat["lam"], cat[star]


def scale_pickles_to_photons(spec_type, mag=0):
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
    
    if type(spec_type) in (list, tuple, np.ndarray) \
        and not hasattr(mag, "__len__"):
        mag = [mag]*len(spec_type)
    else:
        mag = [mag]
    
    Mv = get_stellar_Mv(spec_type)
    if not hasattr(Mv, "__len__"): Mv = [Mv]
    lam, ec = get_pickles_curve(spec_type)
    dlam = (lam[1:] - lam[:-1])
    dlam = np.append(dlam, dlam[-1])
    
    lam *= 1E-4         # convert to um from Ang
    
    # Use Bohlin (2014) to determine the photon flux of a mag 0 A0V star 
    # at 5556 Ang
    F = 3.44E-9 * u.erg / (u.cm**2 * u.s * u.AA)
    E = c.c*c.h/(5556*u.AA)
    ph0 = (F/E).to(1/(u.s * u.cm**2 * u.AA)).value
    
    # 5 Ang/bin * ~962 ph/s * (abs mag + apparent mag)
    ph_factor = []
    for i in range(len(mag)): 
        tmp = dlam * ph0 * 10**(-0.4*(Mv[i] + mag[i]))
        ph_factor += [tmp]    

    # take care of the conversion to ph/s/m2 by multiplying by 1E4
    if type(ec) == (list, tuple):
        for i in len(range(ec)):
            ec[i] *= (lam/0.5556) * ph_factor[i] * 1E4
    else:
        ec *= (lam/0.5556) * ph_factor[0] * 1E4
    
    return lam, ec


def value_at_lambda(lam_i, lam, val):
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
       
    """
    
    i0 = np.where((lam <= lam_i) == True)[0][-1]
    i1 = np.where((lam > lam_i) == True)[0][0]
    
    lam_x = np.array([lam[i0],lam_i,lam[i1]])
    val_i = np.interp(lam_x, lam, val)
    
    return val_i[1]
        

def zero_magnitude_photon_flux(filter_name):
    """
    Return the number of photons for a m=0 star for a certain filter
    
    Parameters
    ----------
    filt : str
        the name of the broadband filter - UBVRIYzJHKKs
        
    Notes
    -----
    units in [ph/s/m2] 
    """
    
    if filter_name not in "UBVRIYzJHKKs": 
        raise ValueError("Filter name must be one of UBVRIYzJHKKs: "+filter_name)
    
    lam, vega = scale_pickles_to_photons("A0V", mag=-0.58)
    
    vraw = ascii.read(os.path.join(__pkg_dir__, "data", 
                                   "TC_filter_"+filter_name+".dat"))
    vlam = vraw[vraw.colnames[0]]
    vval = vraw[vraw.colnames[1]]
    filt = np.interp(lam, vlam, vval)
        
    n_ph = np.sum(vega*filt)
    
    #print("units in [ph/s/m2]")
    return n_ph
    

def SED(spec_type, filter="V", magnitude=0.):
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
    
    if filter not in "UBVRIYzJHKKs": 
        raise ValueError("Filter name must be one of UBVRIYzJHKKs: "+filter)
    
    if type(magnitude) in (list, tuple):
        magnitude = np.asarray(magnitude)
    
    flux_0 = zero_magnitude_photon_flux(filter)
    flux = flux_0 * 10**(-0.4 * magnitude)
        
    lam, star = scale_pickles_to_photons(spec_type)

    vraw = ascii.read(os.path.join(__pkg_dir__, "data", 
                                   "TC_filter_"+filter+".dat"))
    vlam = vraw[vraw.colnames[0]]
    vval = vraw[vraw.colnames[1]]
    filt = np.interp(lam, vlam, vval)
        
    n_ph = np.sum(star*filt)
    
    scale_factor = flux / n_ph
    #print("scale_factor, flux, n_ph, flux_0 [ph/s/m2/bin]")
    #print(scale_factor, flux, n_ph, flux_0, filter)

    return lam,  (scale_factor * star.transpose()).transpose()


def star_grid(n, mag_min, mag_max, filter="K", seperation=1, area=1, spec_type="A0V"):
    """
    Creates a square grid of A0V stars at equal magnitude intervals
    
    Parameters
    ----------
    n : float
        the number of stars in the grid
    mag_min, mag_max : float
        [vega mag] the minimum (brightest) and maximum (faintest) magnitudes for stars in the grid
    filter : str
        broadband filter - B V R I Y z J H K Ks
    seperation : float, optional
        [arcsec] an average speration between the stars in the grid can be specified. Default is 1 arcsec
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
    The weight values are the scaling factors to bring a V=0 A0V spectrum down to the required magnitude for each star
    
    """
    
    if type(mag_min) in (list, tuple, np.ndarray):
        mags = np.asarray(mag_min)
    else:
        if mag_min < mag_max:
            mags = np.linspace(mag_min, mag_max, n)
        elif mag_min < mag_max:
            mags = np.linspace(mag_max, mag_min, n)
        elif mag_min == mag_max:
            mags = np.ones(n) * mag_min
    
    side_len = int(np.sqrt(n)) + (np.sqrt(n) % 1 > 0)
    
    x = seperation * (np.arange(n) % side_len - (side_len - 1) / 2)
    y = seperation * (np.arange(n)// side_len - (side_len - 1) / 2)

    lam, spec = SED(spec_type, filter=filter, magnitude=0)
    if type(spec_type) in (list, tuple):
        ref = np.arange(len(spec_type))
    else:
        ref = np.zeros((n))
    weight = 10**(-0.4*mags) * area
    
    if area != 1:
        units = "ph/s/m2"
    else:
        units = "ph/s"
        
    src = Source(lam=lam, spectra=spec, 
                 x=x, y=y, 
                 ref=ref, weight=weight,
                 units=units)
    
    return src
    
    
def star(mag, filter="K", spec_type="A0V", position=(0,0)):
    """
    Creates a simcado.Source object for a star with a given magnitude
    
    Parameters
    ----------
    mag : float
        magnitude of star
    filt : str
        filter in which the magnitude is given
    spec_type : str, optional
        the spectral type of the star, e.g. "A0V", "G5III"
    
    Returns
    -------
    source : `simcado.Source`
    """
    star = star_grid(1, mag, mag, filter, spec_type=spec_type)
    star.x, star.y = [position[0]], [position[1]]
    return star


def stars(mags, x, y, filter="K", spec_types="A0V"):
    """
    """
    if type(spec_types) in (tuple, list) and len(mags) != len(spec_types):
        raise ValueError("len(mags) != len(spec_types)")
    
    stars = star_grid(len(mags), mags, mags, filter, spec_type=spec_types)
    stars.x, stars.y = x, y
    return stars
    
    
    



