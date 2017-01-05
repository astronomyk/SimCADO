"""
A description of the Chip noise properties and their positions on the Detector

Summary
-------
This module holds three classes: `Detector`, `Chip` and `HXRGNoise`.

`Chip`
Everything to do with photons and electrons happens in the `Chip` class. Each
`Chip` is initialised with a position relative to the centre of the detector
array, a size [in pixels] and a resolution [in arcsec]. Photons fall onto the
`Chip`s and are read out together with the read noise characteristics of the
`Chip`.

`Detector`
The `Detector` holds the information on where the ´Chip`s are placed on the
focal plane. Focal plane coordinates are in [arcsec]. These coordinates are
either read in from a default file or determined by the user. The `Detector`
object is an intermediary - it only passes photons information on the photons
to the `Chip`s. It is mainly a convenience class so that the user can read out
all `Chip`s at the same time.

`HXRGNoise`
This class is borrowed from Berhand Rauscher's script which generates realistic
noise frames for the JWST NIRSpec instrument. NIRSpec uses Hawaii 2RG detectors
but the noise properties scale well up to the H4RG chips that MICADO will use.

Routines
--------
Detector(cmds)
    builds an array of `Chip`s based on a `UserCommands` object
Chip(**kwargs)
    coverts incoming photons into ADUs and adds in read-out noise
HXRGNoise(**kwargs)
    generates realistic detector noise frames

See Also
--------
OpticalTrain, Source

Notes
-----

References
----------
[1] Bernhard Rauscher's HxRG Noise Generator script


Examples
--------
The `Detector` can be used in stand alone mode. In this case it outputs
only the noise that a sealed off detector would generate

```
>>> import simcado
>>> fpa = simcado.Detector(simcado.UserCommands())
>>> fpa.read_out(ouput=True, chips = [0])
...
```

The `Detector` is more useful if we combine it with a `Source` object and an
`OpticalTrain`. Here we create a `Source` object for an open cluster in the LMC
and pass the photons arriving from it through the E-ELT and MICADO. The photons
are then cast onto the detector array. Each `Chip` converts the photons to ADUs
and adds the resulting image to an Astropy `HDUList`. The `HDUList` is then
written to disk.

```
# Create an set of commands, optical train and detector

>>> import simcado
>>> cmds = simcado.UserCommands()
>>> opt_train = simcado.OpticalTrain(cmds)
>>> fpa = simcado.Detector(cmds)

# Pass photons from a 10^4 Msun open cluster in the LMC through to the detector

>>> src = sim.source.source_1E4_Msun_cluster()
>>> src.apply_optical_train(opt_train, fpa)

# Read out the detector array to a FITS file

>>> fpa.read_out(filename="my_raw_image.fits")
```

"""

############################
#       TODO
# - update open, write - remove references to self.params


import os, sys

import warnings, logging
from copy import deepcopy

import multiprocessing as mp

import numpy as np
#from scipy.ndimage.interpolation import zoom

from astropy.io import fits
from astropy.io import ascii as ioascii  # ascii redefines builtin
#from astropy.stats.funcs import median_absolute_deviation as mad

from .utils import __pkg_dir__

from . import spectral as sc
from . import commands
from .nghxrg import HXRGNoise

__all__ = ["Detector", "Chip"]



################################################################################
#                              Detector Objects                                #
################################################################################

class Detector(object):
    """
    Generate a series of `Chip` objects for a focal plane array


    Summary
    -------
    The `Detector` is a holder for the series of `Chip` objects which make up
    the detector array. The main advantage of the `Detector` object is that the
    user can read out all chips in the whole detector array at once. A
    `Detector` is a parameter in the `Source.apply_optical_train()` method.


    Parameters
    ----------
    cmds : UserCommands
        Commands for how to model the Detector
    small_fov : bool, optional
        Default is True. Uses a single 1024x1024 window at the centre of the FoV

    Attributes
    ----------
    cmds : UserCommands
        commands for modelling the detector layout and exposures
    layout : astropy.table.Table
        table of positions and sizes of the chips on the focal plane
    chips : list
        a list of the `Chips` which make up the detector array
    oversample : int
        factor between the internal angular resolution and the pixel FOV
    fpa_res : float
        [mas] field of view of a single pixel
    exptime : float
        [s] exposure time of a single DIT
    tro : float
        [s] time between a single non-destructive readout in up-the-ramp mode
    ndit : int
        number of exposures (DITs)


    Methods
    -------
    read_out(output, filename, chips, **kwargs)
        for reading out the detector array into a FITS file
    open(filename, **kwargs)
        ** not yet implemented **
        Should be moved into a general function for detector.py which returns a
        Detector object after reading in a saved detector file

    write(filename, **kwargs)
        ** not yet implemented **
        Save the Detector object into a FITS file

    Raises
    ------


    See Also
    --------
    Chip, Source, OpticalTrain, UserCommands

    Notes
    -----

    References
    ----------

    Examples
    --------
    Create a `Detector` object
    ```
    >>> import simcado
    >>> my_cmds = simcado.UserCommands()
    >>> my_detector = simcado.Detector(my_cmds)
    ```

    Read out only the first `Chip`
    ```
    >>> my_detector.readout(filename=image.fits, chips=[0])
    ```

    """


    def __init__(self, cmds, small_fov=True):
        # 1. Read in the chip layout
        # 2. Generate chip objects
        # 3. Check if a noise file has been given
            # if not, generate new noise files
            # else: read in the noise file
            # if the noise file has many extensions, choose several random
            #    extensions

        self.cmds = cmds

        if small_fov:
            print("Safety switch is on - Detector(..., small_fov='True')")
            self.layout = ioascii.read(
                """#  id    x_cen    y_cen   x_len   y_len
                   #       arcsec   arcsec   pixel   pixel
                   0        0        0    1024    1024""")
        else:
            self.layout = ioascii.read(self.cmds["FPA_CHIP_LAYOUT"])
        
        self.chips = [Chip(self.layout["x_cen"][i], self.layout["y_cen"][i],
                           self.layout["x_len"][i], self.layout["y_len"][i],
                           self.cmds["SIM_DETECTOR_PIX_SCALE"],
                           self.layout["id"][i])
                      for i in range(len(self.layout["x_cen"]))]

        self.oversample = self.cmds["SIM_OVERSAMPLING"]
        self.fpa_res = self.cmds["SIM_DETECTOR_PIX_SCALE"]
        self.exptime = self.cmds["OBS_EXPTIME"]
        self.ndit    = self.cmds["OBS_NDIT"]
        self.tro     = self.cmds["OBS_NONDESTRUCT_TRO"]
        self._n_ph_atmo   = 0
        self._n_ph_mirror = 0
        self._n_ph_ao = 0
        self.array = None        # defined in method

    def read_out(self, filename=None, to_disk=False, chips=None, **kwargs):
        """
        Simulate the read out process of the detector array

        Summary
        -------
        Based on the parameters set in the `UserCommands` object, the detector
        will read out the images stored on the `Chips` according to the
        specified read out scheme, i.e. Fowler, up-the-ramp, single read, etc.

        Parameters
        ----------
        filename : str
            where the file is to be saved. If `None` the current directory is
            used. Default is `None`
        to_disk : bool
            a flag for where the output should go. If `True` the  `Chip` images
            will be returned to the user (i.e. in an iPython session) as an
            `astropy.fits.HDUList` object. If `False` the `Chip` images will be
            written to a `.fits` file on disk. If no `filename` is specified,
            the output is be called "output.fits". The default is `False`
        chips : int, array-like, optional
            The chip or chips to be read out, based on the detector_layout.dat
            file. Default is the first `Chip` specified in the list, i.e. [0]

        Returns
        -------
        `if to_disk == False:`
            astropy.io.fits.HDUList
        `else:`
            <filename>.fits file

        Keyword Arguments (**kwargs)
        ----------------------------
        **kwargs are used to update the `UserCommands` object which controls
        the `Detector`. Therefore any dictionay keywords can be passed in the
        form of a dictionary, i.e. {"EXPTIME" : 60, "OBS_OUPUT_DIR" : "./"}

        Notes
        -----
        ** In the current release ** the parameter that controls which read
        mode is used is the "SIM_SPEED" command. <3 = full up-the-ramp, 3 <=
        SIM_SPEED < 7 = stacked single read outs. <7 = one single read out for
        the whole NDIT*EXPTIME time span.

        References
        ----------


        Examples
        --------
        """

        #removed kwargs
        self.cmds.update(kwargs)

        if filename is not None:
            to_disk = True

        if filename is None and to_disk is True:
            if self.cmds["OBS_OUTPUT_DIR"] is None:
                self.cmds["OBS_OUTPUT_DIR"] = "./output.fits"
            filename = self.cmds["OBS_OUTPUT_DIR"]

        if chips is not None:
            if not hasattr(chips, "__len__"):
                ro_chips = [chips]
            else:
                ro_chips = chips
        elif chips is None:
            ro_chips = np.arange(len(self.chips))
        else:
            raise ValueError("Something wrong with `chips`")

        hdulist = fits.HDUList()
        if len(ro_chips) > 1:
            primary_hdu = fits.PrimaryHDU()
            for key in self.cmds.cmds:
                val = self.cmds.cmds[key]
                
                if isinstance(val, (sc.TransmissionCurve, sc.EmissionCurve, 
                                    sc.UnityCurve, sc.BlackbodyCurve)):
                    val = val.params["filename"]
                
                if isinstance(val, str) and len(val) > 35:
                    val = "... " + val[-35:]
                
                try:
                    primary_hdu.header["HIERARCH "+key] = val
                except NameError:   # any other exceptions possible?
                    pass
            hdulist.append(primary_hdu)

        for i in ro_chips:
            ######
            # Put in a catch here so that only the chips specified in "chips"
            # are read out
            ######
            print("Reading out chip", self.chips[i].id)
            array = self.chips[i].read_out(self.cmds)

            ## TODO: transpose is just a hack - need to make sure
            ##       x and y are handled correctly throughout SimCADO
            thishdu = fits.ImageHDU(array.T)

            thishdu.header["EXTNAME"] = ("CHIP_{:02d}".format(self.chips[i].id),
                                         "Chip ID")

            thishdu.header["CTYPE1"] = "LINEAR"
            thishdu.header["CUNIT1"] = "arcsec"
            thishdu.header["CRVAL1"] = 0.
            thishdu.header["CRPIX1"] = (self.chips[i].naxis1 //2 -
                                        self.chips[i].x_cen / self.chips[i].pix_res)
            thishdu.header["CDELT1"] = (self.chips[i].pix_res,
                                        "[arcsec] Pixel scale")

            # axis 2
            thishdu.header["CTYPE2"] = "LINEAR"
            thishdu.header["CUNIT2"] = "arcsec"
            thishdu.header["CRVAL2"] = 0.
            thishdu.header["CRPIX2"] = (self.chips[i].naxis2 //2 -
                                        self.chips[i].y_cen / self.chips[i].pix_res)
            thishdu.header["CDELT2"] = (self.chips[i].pix_res,
                                        "[arcsec] Pixel resolution")

            # possible rotation
            thishdu.header["PC1_1"] = 1.
            thishdu.header["PC1_2"] = 0.
            thishdu.header["PC2_1"] = 0.
            thishdu.header["PC2_2"] = 1.

            thishdu.header["CHIP_ID"] = (self.chips[i].id, "Chip ID")
            thishdu.header["BUNIT"] = ("ph/s", "")
            thishdu.header["EXPTIME"] = (self.exptime, "[s] Exposure time")
            thishdu.header["NDIT"] = (self.ndit, "Number of exposures")
            thishdu.header["TRO"] = (self.tro,
                                     "[s] Time between non-destructive readouts")

            for key in self.cmds.cmds:
                val = self.cmds.cmds[key]
                if isinstance(val, str): 
                    if len(val) > 35:
                        val = "... " + val[-35:]
                try:
                    thishdu.header["HIERARCH "+key] = val
                except NameError:   # any other exceptions possible?
                    pass
                except ValueError:
                    warnings.warn("ValueError - Couldn't add keyword: "+key)
            hdulist.append(thishdu)
            # hdulist = fits.HDUList(hdus)

        if not to_disk:
            return hdulist
        else:
            hdulist.writeto(filename, clobber=True, checksum=True)
    
    
    def write(self, filename=None, **kwargs):
        """
        Write a `Detector` object out to a FITS file


        Summary
        -------
        Writes the important information containeed in a `Detector` object into
        FITS file for later use. The main information written out include: the
        layout of the detector chips, any pixel maps associated with the
        detector chips, a linearity curve and a QE curve for the chips.


        Parameters
        ----------
        filename : str, optional
            path to the FITS file where the `Detector` object is stored. If
            `filename=None` (by default), the file written is `./detector.fits`


        Returns
        -------
        None

        Keyword Arguments (**kwargs)
        ----------------------------

        Examples
        --------
        """

        self.params.update(kwargs)
        if filename is None:
            filename = self.params["OBS_OUTPUT_DIR"]
        if filename is None:
            raise ValueError("No output path was specified. " + \
                             "Use either filename= or HXRG_OUTPUT_PATH=")

        hdu = fits.PrimaryHDU(self.array)
        hdu.header["PIX_RES"] = self.pix_res
        hdu.header["EXPTIME"] = self.exptime
        hdu.header["GAIN"]    = self.params["FPA_GAIN"]
        hdu.header["SIMCADO"] = "FPA_NOISE"

        try:
            hdu.writeto(filename, clobber=True, checksum=True)
        except OSError:
            warnings.warn(filename+" exists and is busy. OS won't let me write")


def open(self, filename):
    """
    Opens a saved `Detector` file.


    Summary
    -------
    ** Not yet implemented **
    ** Should be moved outside of `Detector` and called with
    `detector.open()` **

    Detector objects can be saved to FITS file and read back in for later
    simulations.

    Parameters
    ----------
    filename : str
        path to the FITS file where the `Detector` object is stored

    Returns
    -------
    `simcado.Detector` object

    Examples
    --------
    """

    if not os.path.exists(filename):
        raise FileNotFoundError(filename + " doesn't exist")

    with fits.open(filename) as fp1:
        self.params.update(fp1[0].header)
        self.array = fp1[0].data

            
            

def plot_detector_layout(detector):
    """Plot the detector layout. NOT FINISHED """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ValueError("matplotlib can't be found")

    plt.figure(figsize=(10, 10))
    clr = ["g"]

    for i in range(len(detector.chips)):
        chip = detector.chips[i]
        plt.plot((chip.x_min, chip.x_max), (chip.y_min, chip.y_min),
                 c=clr[i%len(clr)])
        plt.plot((chip.x_min, chip.x_max), (chip.y_max, chip.y_max),
                 c=clr[i%len(clr)])
        plt.plot((chip.x_min, chip.x_min), (chip.y_min, chip.y_max),
                 c=clr[i%len(clr)])
        plt.plot((chip.x_max, chip.x_max), (chip.y_min, chip.y_max),
                 c=clr[i%len(clr)])
        plt.text(chip.x_cen, chip.y_cen, chip.id, fontsize=14)
        plt.xlabel("Distance [arcsec]", fontsize=14)
        plt.ylabel("Distance [arcsec]", fontsize=14)

    plt.show()


def plot_detector(detector):
    """
    Plot the contents of a detector array

    Parameters
    ----------
    detector : simcado.Detector
        The detector object to be shown
    """

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    plt.figure(figsize=(15*2, 13.75*2))

    x0 = np.min([i.x_min for i in detector.chips])
    x1 = np.max([i.x_max for i in detector.chips])
    y0 = np.min([i.y_min for i in detector.chips])
    y1 = np.max([i.y_max for i in detector.chips])
    w = (detector.chips[0].x_max - detector.chips[0].x_min)/(x1 - x0)
    h = (detector.chips[0].y_max - detector.chips[0].y_min)/(y1 - y0)

    for chip in detector.chips:
        s = plt.axes([(chip.x_min - x0)/(x1 - x0), (chip.y_min - y0)/(y1 - y0),
                      w, h])
        s.set_xticks([])
        s.set_yticks([])
        if chip.array is not None:
            s.imshow(np.rot90(chip.array - np.min(chip.array)), norm=LogNorm(),
                     cmap="Greys", vmin=1)

    plt.show()


################################################################################
#                              Chip Objects                                    #
################################################################################


class Chip(object):
    """
    Holds the "image" as seen my a single chip in the focal plane


    Summary
    -------
    The `Chip` object contains information on where it is located in the focal
    plane array. The method `<Source>.apply_optical_train()` passes an image of
    the on-sky object to each `Chip`. THis image is resampled to the `Chip`
    pixel scale. Each `Chip` holds the "ideal" image as an array of expectraion
    values for the level of photons arriving during an EXPTIME. The `Chip` then
    adds detector noise and other characteristics to the image when
    <Detector>.readout() is called.


    Parameters
    ----------
    x_cen, y_cen : float
        [arcsec] the coordinates of the centre of the chip relative to the
        centre of the focal plane
    x_len, y_len : int
        the number of pixels per dimension
    pix_res : float
        [arcsec] the field of view per pixel
    id : int
        an indetification number for the chip (assuming they are not correctly
        ordered)

    Attributes
    ----------
    x_cen, y_cen : float
        [arcsec] the coordinates of the centre of the chip relative to the
        centre of the focal plane
    nasxis1, naxis2 : int
        the number of pixels per dimension
    pix_res : float
        [arcsec] the field of view per pixel
    chipid : int, optional
        the id of the chip relative to the others on the detector array. Default is `None`
    dx, dy : float
        [arcsec] half of the field of view of each chip
    x_min, x_max, y_min, y_max : float
        [arcsec] the borders of the chip realtive to the centre of the focal plane
    array : np.ndarray
        an array for holding the signal registered by the `Chip`


    Methods
    -------
    add_signal(signal)
        adds signal to `.array`. The signal should be the same dimensions as
        `Chip.array`
    add_uniform_background(emission, lam_min, lam_max, output=False)
        adds a constant to the signal in `.array`. The background level is found
        by integrating the the `emission` curve between `lam_min` and `lam_max`.
        It output is set to `True`, an image with the same dimensions as
        `.array` scaled to the bacground flux is returned
    apply_pixel_map(pixel_map_path=None, dead_pix=None, max_well_depth=1E5)
        applies a mask to `.array` representing the position of the current
        "hot" and "dead" pixels / lines
    reset()
        resets the signal on the `Chip` to zero. In future releases, an
        imlementation of the persistence characteristics of the detector will
        go here.


    Raises
    ------

    See Also
    --------
    Detector, Source, UserCommands, OpticalTrain

    Examples
    --------
    """

    def __init__(self, x_cen, y_cen, x_len, y_len, pix_res, chipid=None):

        self.x_cen  = x_cen
        self.y_cen  = y_cen
        self.naxis1 = x_len
        self.naxis2 = y_len
        self.pix_res = pix_res
        self.id     = chipid      # id is built-in, should not be redefined

        dx = (x_len // 2) * pix_res
        dy = (y_len // 2) * pix_res
        self.x_min = x_cen - dx
        self.x_max = x_cen + dx
        self.y_min = y_cen - dy
        self.y_max = y_cen + dy

        self.array = None


    def add_signal(self, signal):
        """
        Add a 2D array of photon signal to the Chip

        Summary
        -------
        Add some signal photons to the detector array. Input units are expected
        to be [ph/s/pixel]

        Parameters
        ----------
        signal : np.ndarray
            [ph/pixel/s] photon signal. `signal` should have the same dimensions
            as the `array`

        Returns
        -------
        None

        Raises
        ------

        Examples
        --------

        """
        if isinstance(signal, np.ndarray):
            if signal.shape[0] == self.naxis1 and \
               signal.shape[1] == self.naxis2:
                if self.array is None:
                    self.array = signal
                else:
                    self.array += signal
            else:
                raise ValueError(str(signal.shape) + " != " + \
                                 str((self.naxis1, self.naxis2)))
        elif not hasattr(signal, "__len__"):
            self.array += signal


    def add_uniform_background(self, emission, lam_min, lam_max, output=False):
        """
        Add a uniform background

        Summary
        -------
        Take an EmissionCurve and some wavelength boundaries, lam_min lam_max,
        and sum up the photons in between. Add those to the source array.

        Parameters
        ----------
        - emission_curve: EmissionCurve object with background emission photons
        - lam_min, lam_max: the wavelength limits

        Optional keywords:
        - output: [False, True] if output is True, the BG emission array is
                  returned

        Output is in [ph/s/pixel]
        """

        if isinstance(emission, sc.EmissionCurve):
            bg_photons = emission.photons_in_range(lam_min, lam_max)
        elif isinstance(emission, (float, int)):
            bg_photons = emission
        else:
            bg_photons = 0
            warnings.warn("type(emission) invalid. No background added")

        if output is True:
            return bg_photons * np.ones(self.array.shape, dtype=np.float32)
        else:
            self.array += bg_photons


    def apply_pixel_map(self, pixel_map_path=None, dead_pix=None,
                        max_well_depth=1E5):
        """
        adds "hot" and "dead" pixels to the array


        Summary
        -------
        applies a mask to `.array` representing the positions of the current
        "hot" and "dead" pixels / lines. The method either reads in a FITS file
        with locations of these pixels, or generates a series of random
        coordinates and random weights for the pixels.

        Parameters
        ----------
        pixel_map_path : str
            path to the FITS file. Default is None
        dead_pix : int
            [%] the percentage of dead or hot pixels on the chip - only used if
            pixel_map_path = None. Default is `None`.
        max_well_depth : 1E5


        Returns
        -------
        None

        Raises
        ------

        See Also

        Examples
        --------
        """

        try:
            pixel_map = fits.getdata(pixel_map_path)
            if self.array.shape != pixel_map.shape:
                raise ValueError("pixel_map.shape != detector_array.shape")
            self.array += pixel_map * max_well_depth
        except ValueError:
            if dead_pix is not None:
                n = int(self.naxis1 * self.naxis2 * dead_pix / 100)
                x = np.random.randint(self.naxis1, size=n)
                y = np.random.randint(self.naxis2, size=n)
                z = np.random.random(n)
                self.array[x, y] += z * max_well_depth
            else:
                raise ValueError("Couldn't apply pixel_map")


    def _apply_saturation(self, arr):
        """
        Cap all pixels that are above the well depth.

        Parameters
        ----------
        x  :  type [, optional [, {set values} ]]
            Description of `x`. [(Default value)]

        Returns
        -------

        Examples
        --------
        """
        # TODO: apply a linearity curve and shift excess light into !!
        # neighbouring pixels !!

        max_val = self.params["FPA_WELL_DEPTH"]
        arr[arr > max_val] = max_val
        return arr


    def reset(self):
        """
        <One-line summary goes here>

        Summary
        -------

        Parameters
        ----------
        x  :  type [, optional [, {set values} ]]
            Description of `x`. [(Default value)]

        Returns
        -------
        """

        ### TODO - add in persistence
        self.array = None


    def read_out(self, cmds):
        """
        Readout the detector array

        Summary
        -------

        Parameters
        ----------
        cmds : simcado.UserCommands
            Commands for how to read out the chip

        Returns
        -------
        out_array : np.ndarray
            image of the chip read out

        Raises
        ------

        Examples
        --------
        """

        ###############################################
        #!!!!!!!!!!! TODO - add dark strom !!!!!!!!!!!#

        dit      = cmds["OBS_EXPTIME"]
        ndit     = int(cmds["OBS_NDIT"])
        #tro      = cmds["OBS_NONDESTRUCT_TRO"]
        #max_byte = cmds["SIM_MAX_RAM_CHUNK_GB"] * 2**30
        dark     = cmds["FPA_DARK_MEDIAN"]

        if self.array is None:
            self.array = np.zeros((self.naxis1, self.naxis2), dtype=np.float32)

        # At this point, the only negatives come from the convolution.
        # Remove them for the Poisson process
        self.array[self.array < 0] = 0

        out_array = np.zeros(self.array.shape, dtype=np.float32)

        ############## TO DO add in the different read out modes ###############
        # if cmds["SIM_SPEED"] <= 3:
            # for n in range(ndit):
                # out_array += self._read_out_uptheramp(self.array, dit, tro,
                                                   # max_byte)
                # out_array += dark * dit

        # elif cmds["SIM_SPEED"] > 3 and cmds["SIM_SPEED"] <= 7:
            # for n in range(ndit):
                # out_array += self._read_out_fast(self.array, dit)

            # out_array += dark * dit * ndit

        #elif cmds["SIM_SPEED"] > 7:
        #######################################################################

        out_array = self._read_out_superfast(self.array, dit, ndit)
        
        # apply the linearity curve
        if cmds["FPA_LINEARITY_CURVE"] is not None:
            fname = cmds["FPA_LINEARITY_CURVE"]

            data = ioascii.read(fname)
            real_cts = data[data.colnames[0]]
            measured_cts =  data[data.colnames[1]]

            out_array = np.interp(out_array.flatten(), 
                              real_cts, measured_cts).reshape(out_array.shape)
            out_array = out_array.astype(np.float32)
        
        #### TODO #########
        # add read out noise for every readout
        # add a for loop to _read_noise where n random noise frames are added
        # based on the size of the noise cube
        ro = self._read_noise_frame(cmds)
        for i in range(1, ndit):
            ro += self._read_noise_frame(cmds) + dark * dit

        out_array += ro

        if cmds["OBS_REMOVE_CONST_BG"].lower() == "yes":
            min_val = np.min(out_array)
            out_array -= min_val

        return out_array

    ## TODO: What to do if dit = min_dit (single read)?
    ## TODO: Make breaking up into memory chunks more flexible?
    def _read_out_uptheramp(self, image, dit, tro=1.3, max_byte=2**30):
        """
        Test readout onto a detector using cube model

        Parameters
        ==========
        image :
            a 2D image to be mapped onto the detector. Units are [ph/s/pixel]
        dit :
            integration time [s]
        tro :
            time for a single non-destructive read (default: 1.3 seconds)

        Optional Parameters
        ===================
        max_byte :
            the largest possible chunk of memory that can be used for computing
            the sampling slope

        This function builds an intermediate cube of dimensions (nx, ny, nro) with a
        layer for  each non-destructive read.

        Output is given in [ph/pixel]

        """
        nx, ny = image.shape

        nro = np.int(dit / tro)
        tpts = (1 + np.arange(nro)) * tro

        img_byte = image.nbytes
        pix_byte = img_byte / (nx * ny)

        #max_byte =            ## TODO: arbitrary, function parameter?
        max_pix = max_byte / pix_byte

        #cube_megabyte = img_byte * nro / 2**20
        #print("Full cube  has {0:.1f} Megabytes".format(cube_megabyte))

        ny_cut = np.int(max_pix / (nx * nro))
        if ny_cut >= ny:
            ny_cut = ny
        #print("Cut image to ny={0:d} rows".format(ny_cut))

        slope = np.zeros(image.shape)
        cube = np.zeros((nx, ny_cut, nro))

        y1 = 0
        while y1 < ny:

            y2 = y1 + ny_cut
            if y2 > ny:
                y2 = ny
                ny_cut = ny - y1
                try:
                    del cube
                    cube = np.zeros((nx, ny_cut, nro))
                except NameError:
                    pass

            ## Fill the cube with Poisson realization, individual reads
            for i in range(nro):
                cube[:, :, i] = np.random.poisson(image[:, y1:y2] * tro)

            ## Build the ramp
            sumcube = cube.cumsum(axis=2)

            ## determine the slope using explicit formula calculated over cube
            Sx = tpts.sum()
            Sxx = (tpts * tpts).sum()
            Sy = np.sum(sumcube, axis=2)
            Sxy = np.sum(sumcube * tpts, axis=2)

            slope[:, y1:y2] = (nro * Sxy - Sx * Sy) / (nro * Sxx - Sx * Sx)

            ## Move to next slice
            y1 = y2

        # return values are [ph/pixel]
        return slope * dit


    def _read_out_poisson(self, image, dit, ndit):
        """
        <One-line summary goes here>

        Parameters
        ----------
        x : type [, optional [, {set values} ]]
            Description of `x`. [(Default value)]

        Returns
        -------

        Examples
        --------
        """
        image2 = image * dit
        image2[image2 > 2.14E9] = 2.14E9

        im_st = np.zeros(np.shape(im))
        for i in range(ndit):
            im_st += np.random.poisson(image2)
            
        return im_st.astype(np.float32)

        
    def _read_out_fast(self, image, dit):
        """
        <One-line summary goes here>

        Parameters
        ----------
        x  :  type [, optional [, {set values} ]]
            Description of `x`. [(Default value)]

        Returns
        -------

        Examples
        --------
        """

        image2 = image * dit
        image2[image2 > 2.14E9] = 2.14E9
        
        return np.random.poisson(image2 * dit)


    def _read_out_superfast(self, image, dit, ndit):
        """
        <One-line summary goes here>

        Summary
        -------

        Parameters
        ----------
        x  :  type [, optional [, {set values} ]]
            Description of `x`. [(Default value)]

        Returns
        -------

        Examples
        --------
        """

        exptime = dit * ndit
        image2 = image * exptime
        image2[image2 > 2.14E9] = 2.14E9
        return np.random.poisson(image2).astype(np.float32)


    def _read_noise_frame(self, cmds):
        """
        Read in read-out-noise from the FITS file specified by FPA_NOISE_PATH

        Parameters
        ----------
        x  :  type [, optional [, {set values} ]]
            Description of `x`. [(Default value)]

        Returns
        -------

        Raises
        ------

        Examples
        --------
        """

        if cmds["FPA_USE_NOISE"].lower() == "no":
            return np.zeros((self.naxis1, self.naxis2))

        if "gen" in cmds["FPA_NOISE_PATH"].lower():
            if cmds["HXRG_OUTPUT_PATH"] is not None:
                _generate_hxrg_noise(self.naxis1, self.naxis2, cmds)
                tmp = fits.getdata(cmds["HXRG_OUTPUT_PATH"])
                return tmp[:self.naxis1, :self.naxis2]
            else:
                return _generate_hxrg_noise(self.naxis1, self.naxis2, cmds)

        elif cmds["FPA_NOISE_PATH"] is not None:
            n = len(fits.info(cmds["FPA_NOISE_PATH"], False))
            layer = np.random.randint(n)
            tmp = fits.getdata(cmds["FPA_NOISE_PATH"], layer)
            return tmp[:self.naxis1, :self.naxis2]
        else:
            return np.zeros((self.naxis1, self.naxis2))


    def __array__(self):
        return self.array

    def __mul__(self, x):
        psf_new = deepcopy(self)
        return psf_new.array * x

    def __add__(self, x):
        psf_new = deepcopy(self)
        return psf_new.array + x

    def __sub__(self, x):
        psf_new = deepcopy(self)
        return psf_new.array - x

    def __rmul__(self, x):
        return self.__mul__(x)

    def __radd__(self, x):
        return self.__add__(x)

    def __rsub__(self, x):
        psf_new = deepcopy(self)
        return x - psf_new.array

    def __imul__(self, x):
        return self.__mul__(x)

    def __iadd__(self, x):
        return self.__add__(x)

    def __isub__(self, x):
        return self.__sub__(x)


def _generate_hxrg_noise(cmds):
#def _generate_hxrg_noise(naxis1, naxis2, cmds):
    """
    Generate a read noise frame using a UserCommands object

    Create a detector noise array using Bernard Rauscher's NGHxRG tool

    Parameters
    ----------
    cmds : simcado.UserCommands

    Returns
    -------

    Examples
    --------
    """
    import multiprocessing as mp
    
    #if len(kwargs) > 0 and self.verbose: print("updating ",kwargs)
    #self.params.update(kwargs)
    print("Generating a new chip noise array")
    print(mp.current_process())
    # HXRG needs a pca file to run. Work out what a PCA file means!!
    ng_h4rg = HXRGNoise(naxis1=4096,
                        naxis2=4096,
                        naxis3=1,
                        n_out=cmds["HXRG_NUM_OUTPUTS"],
                        nroh=cmds["HXRG_NUM_ROW_OH"],
                        pca0_file=cmds["HXRG_PCA0_FILENAME"],
                        verbose=cmds["SIM_VERBOSE"])

    # Make a noise file
    noise = ng_h4rg.mknoise(o_file=cmds["HXRG_OUTPUT_PATH"],
                            rd_noise=cmds["FPA_READOUT_MEDIAN"],
                            pedestal=cmds["HXRG_PEDESTAL"],
                            c_pink=cmds["HXRG_CORR_PINK"],
                            u_pink=cmds["HXRG_UNCORR_PINK"],
                            acn=cmds["HXRG_ALT_COL_NOISE"])

    return noise


def make_noise_cube(num_layers=25, filename="FPA_noise.fits", multicore=True):
    """
    Create a large noise cube with many separate readout frames.

    Note:
    Each frame take about 15 seconds to be generated. The default value of
    25 frames will take around six minutes depending on your computer's
    architecture.

    Parameters
    ----------
    num_layers : int, optional
        the number of separate readout frames to be generated. Default is 25
    filename : str, optional
        The filename for the FITS cube. Default is "FPA_noise.fits"
    multicore : bool, optional
        If you're not using windows, this allows the process to use all
        available cores on your machine to speed up the process. Default is True

    Notes
    -----
    multicore doesn't work - fix it

    """
    
    if sys.version_info.major >= 3:
        print("Sorry, but this only works in Python 3 and above. \
           See the SimCADO FAQs for work-around options")
        return None
    
    
    cmds = commands.UserCommands()
    cmds["FPA_NOISE_PATH"] = "generate"
    cmds["FPA_CHIP_LAYOUT"] = "default"

    layout = ioascii.read(cmds.cmds["FPA_CHIP_LAYOUT"])
    naxis1, naxis2 = layout["x_len"][0], layout["y_len"][0]

    if "Windows" in os.environ.get('OS',''):
        multicore = False

    if __name__ == "__main__" and multicore:
        pool = mp.Pool(processes=mp.cpu_count()-1)
        frames = pool.map(_generate_hxrg_noise, (cmds)*num_layers)
    else:
        frames = [_generate_hxrg_noise(cmds) \
                  for i in range(num_layers)]
    
    hdu = fits.HDUList([fits.PrimaryHDU(frames[0])] + \
        [fits.ImageHDU(frames[i]) \
        for i in range(1, num_layers)])
    
    if filename is None:
        return hdu
    else:
        hdu.writeto(filename, clobber=True, checksum=True)


def install_noise_cube(n=9):
    """
    Install a noise cube in the package directory

    Parameters
    ----------
    n : int, optional
        number of layers.

    Warning
    -------
    Each layer is ~64MB, default is 9 layers (~600MB). If you have less than
    1 GB on the drive where your Python installation is. Be careful!
    """

    if sys.version_info.major >= 3:
        print("WARNING - this process can take up to 10 minutes. Fear not!")
        hdu = make_noise_cube(n, filename=None)
        filename = os.path.join(__pkg_dir__, "data", "FPA_noise.fits")
        hdu.writeto(filename, clobber=True, checksum=True)
        print("Saved noise cube with", n, "layers to the package directory:")
        print(filename)
    else:
        print("Sorry, but this only works in Python 3 and above. \
               See the SimCADO FAQs for work-around options")
    
    