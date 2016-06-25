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

>>> src = sim.optics_utils.source_1E4_Msun_cluster()
>>> src.apply_optical_train(opt_train, fpa)

# Read out the detector array to a FITS file

>>> fpa.read_out(filename="my_raw_image.fits")
```

"""

############################
#       TODO
# - update open, write - remove references to self.params


import os, warnings, datetime

import numpy as np
from scipy.ndimage.interpolation import zoom

from astropy.io import fits, ascii

try:
    import simcado.spectral as sc
except:
    import spectral as sc

__all__ = ["Detector", "Chip"]


    
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


    def __init__(self, cmds, small_fov=False):
        # 1. Read in the chip layout
        # 2. Generate chip objects
        # 3. Check if a noise file has been given
            # if not, generate new noise files
            # else: read in the noise file
            # if the noise file has many extensions, choose several random 
            #    extensions
        
        self.cmds = cmds
        
        if small_fov:
            self.layout = ascii.read("""#  id    x_cen    y_cen   x_len   y_len
                                        #       arcsec   arcsec   pixel   pixel
                                         0        0        0    1024    1024""")
        else:
            self.layout = ascii.read(self.cmds["FPA_CHIP_LAYOUT"])
        self.chips  = [Chip(self.layout["x_cen"][i], self.layout["y_cen"][i], 
                            self.layout["x_len"][i], self.layout["y_len"][i], 
                            self.cmds["SIM_DETECTOR_PIX_SCALE"], 
                            self.layout["id"][i]) 
                       for i in range(len(self.layout["x_cen"]))]
        
        self.oversample = self.cmds["SIM_OVERSAMPLING"]
        self.fpa_res = self.cmds["SIM_DETECTOR_PIX_SCALE"]
        self.exptime = self.cmds["OBS_EXPTIME"]
        self.ndit    = self.cmds["OBS_NDIT"]
        self.tro     = self.cmds["OBS_NONDESTRUCT_TRO"]
            

    def read_out(self, filename=None, to_disk=False, chips=None):
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
        `if output == True:`
            astropy.io.fits.HDUList
        `else:`
            <filename>.fits file 
        


        Keyword Arguments (**kwargs)
        ----------------------------
        **kwargs are used to update the `UserCommands` object which controls 
        the `Detector`. Therefore any dictionay keywords can be passed in the 
        form of a dictionary, i.e. {"EXPTIME" : 60, "OBS_OUPUT_DIR" : "./"}

        Raises
        ------


        See Also
        --------
        

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
        #self.cmds.update(kwargs)
        
        if filename is None and to_disk is True: 
            if self.cmds["OBS_OUTPUT_DIR"] is None:
                self.cmds["OBS_OUTPUT_DIR"] = "./output.fits"
            filename = self.cmds["OBS_OUTPUT_DIR"]
        
        if chips is not None and not hasattr(chips, "__len__"):
            ro_chips = [chips]
        elif chips is not None and hasattr(chips, "__len__"):
            ro_chips = chips
        elif chips is None:
            ro_chips = np.arange(len(self.chips))
        else:
            raise ValueError("Something wrong with `chips`")
        
        pri_hdu_flag = False
        for i in ro_chips:
            ######
            # Put in a catch here so that only the chips specified in "chips"
            # are read out
            ######
            
            
            array = self.chips[i].read_out(self.cmds)
                        
            if pri_hdu_flag == False:
                hdus = [fits.PrimaryHDU(array)]
                pri_hdu_flag = True
            else:
                hdus += [fits.ImageHDU(array)]
            
            hdus[i].header["CDELT1"] = (self.chips[i].pix_res, 
                                                    "[arcsec] Pixel resolution")
            hdus[i].header["CDELT2"] = (self.chips[i].pix_res,
                                                    "[arcsec] Pixel resolution")
            hdus[i].header["CRVAL1"] = (self.chips[i].x_cen, 
                        "[arcsec] central pixel relative to detector centre")
            hdus[i].header["CRVAL2"] = (self.chips[i].y_cen, 
                        "[arcsec] central pixel relative to detector centre")
            hdus[i].header["CRPIX1"] = (self.chips[i].naxis1 // 2, 
                                                                "central pixel")
            hdus[i].header["CRPIX2"] = (self.chips[i].naxis2 // 2, 
                                                                "central pixel")
            hdus[i].header["CHIP_ID"] = (self.chips[i].id, "Chip ID")

            hdus[i].header["BUNIT"] = ("ph/s", "")
            hdus[i].header["EXPTIME"] = (self.exptime, "[s] Exposure time")
            hdus[i].header["NDIT"]  = (self.ndit, "Number of exposures")
            hdus[i].header["TRO"]   = (self.tro, 
                                    "[s] Time between non-destructive readouts")

        hdulist = fits.HDUList(hdus)
        
        if to_disk is False:
            return hdulist
        else:
            hdulist.writeto(filename, clobber=True)
        
            
    def open(self, filename, **kwargs):
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


        Keyword Arguments (**kwargs)
        ----------------------------


        Raises
        ------


        See Also
        --------


        Notes
        -----


        References
        ----------


        Examples
        --------
        """
     
        if not os.path.exists(filename):
            raise ValueError(filename + " doesn't exist")

        f = fits.open(filename)
        self.params.update(f[0].header)
        self.array = f[0].data
        f.close()

        
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

        Raises
        ------

        See Also
        --------

        Notes
        -----

        References
        ----------

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
        hdu.header["SIMCADO"]= "FPA_NOISE"

        try:
            hdu.writeto(filename, clobber=True)
        except:
            warnings.warn(filename+" exists and is busy. OS won't let me write")

    
def plot_detector_layout(detector):
    """Plot the detector layout. NOT FINISHED """
    try:
        import matplotlib.pyplot as plt
    except:
        raise ValueError("matplotlib can't be found")
    
    plt.figure(figsize=(10,10))
    clr = ["g"]

    for i in range(9):
        chip = detector.chips[i]
        plt.plot((chip.x_min,chip.x_max), (chip.y_min,chip.y_min),c=clr[i%len(clr)])
        plt.plot((chip.x_min,chip.x_max), (chip.y_max,chip.y_max),c=clr[i%len(clr)])
        plt.plot((chip.x_min,chip.x_min), (chip.y_min,chip.y_max),c=clr[i%len(clr)])
        plt.plot((chip.x_max,chip.x_max), (chip.y_min,chip.y_max),c=clr[i%len(clr)])
        plt.text(chip.x_cen,chip.y_cen,str(i),fontsize=14)
        plt.xlabel("Distance [arcsec]", fontsize=14); plt.ylabel("Distance [arcsec]", fontsize=14)

    plt.show()

        
    
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
    id : int, optional
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
    reset_chip()
        resets the signal on the `Chip` to zero. In future releases, an 
        imlementation of the persistence characteristics of the detector will 
        go here.
       

    Raises
    ------

    See Also
    --------
    Detector, Source, UserCommands, OpticalTrain

    Notes
    -----


    References
    ----------


    Examples
    --------
    """

    def __init__(self, x_cen, y_cen, x_len, y_len, pix_res, id=None):
        
        self.x_cen  = x_cen
        self.y_cen  = y_cen
        self.naxis1 = x_len
        self.naxis2 = y_len
        self.pix_res = pix_res
        self.id     = id
        
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

        See Also
        --------

        Notes
        -----

        References
        ----------

        Examples
        --------
        
        """
        if type(signal) == np.ndarray:
            if signal.shape[0] == self.naxis1 and signal.shape[1] == self.naxis2:
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
        <One-line summary goes here>


        Summary
        -------


        Parameters
        ----------
        x : type [, optional [, {set values} ]]
            Description of `x`. [(Default value)]

        Returns
        -------


        Keyword Arguments (**kwargs)
        ----------------------------


        Raises
        ------


        See Also
        --------


        Notes
        -----


        References
        ----------


        Examples
        --------
        """

        """
        Take an EmissionCurve and some wavelength boundaries, lam_min lam_max,
        and sum up the photons in between. Add those to the source array.

        Keywords:
        - emission_curve: EmissionCurve object with background emission photons
        - lam_min, lam_max: the wavelength limits

        Optional keywords:
        - output: [False, True] if output is True, the BG emission array is
                  returned
                  
        Output is in [ph/s/pixel]
        """
        
        if type(emission) == sc.EmissionCurve:
            bg_photons = emission.photons_in_range(lam_min, lam_max)        
        elif type(emission) in (float, int):
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
        --------

        Notes
        -----

        References
        ----------

        Examples
        --------
        """
 
        try:
            pixel_map = fits.getdata(pixel_map_path)
            if self.array.shape != pixel_map.shape:
                raise ValueError("pixel_map.shape != detector_array.shape")
            self.array += pixel_map * max_well_depth
        except:
            if dead_pix is not None:
                n = int(self.naxis1 * self.naxis2 * dead_pix / 100)
                x = np.random.randint(self.naxis1, size=n)
                y = np.random.randint(self.naxis2, size=n)
                z = np.random.random(n)
                self.array[x,y] += z * max_well_depth
            else:
                raise ValueError("Couldn't apply pixel_map")


    def _apply_saturation(self, arr):
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


        Keyword Arguments (**kwargs)
        ----------------------------


        Raises
        ------


        See Also
        --------


        Notes
        -----


        References
        ----------


        Examples
        --------
        """

        """
        Cap all pixels that are above the well depth.
        !! TODO: apply a linearity curve and shift excess light into !!
        !! neighbouring pixels !!
        """
        max_val = self.params["FPA_WELL_DEPTH"]
        arr[arr > max_val] = max_val
        return arr
                                 
                                 
    def reset_chip(self):
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


        Keyword Arguments (**kwargs)
        ----------------------------


        Raises
        ------


        See Also
        --------


        Notes
        -----


        References
        ----------


        Examples
        --------
        """

        """
        """
        ### TODO - add in persistence 
        self.array = None


    def read_out(self, cmds):
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


        Keyword Arguments (**kwargs)
        ----------------------------


        Raises
        ------


        See Also
        --------


        Notes
        -----


        References
        ----------


        Examples
        --------
        """

        """
        Readout the detector array
        """

        ###############################################
        #!!!!!!!!!!! TODO - add dark strom !!!!!!!!!!!#
        
        dit     = cmds["OBS_EXPTIME"]
        ndit    = int(cmds["OBS_NDIT"])
        tro     = cmds["OBS_NONDESTRUCT_TRO"]
        max_byte = cmds["SIM_MAX_RAM_CHUNK_GB"] * 2**30
        dark    = cmds["FPA_DARK_MEDIAN"]
        
        if self.array is None:
            self.array = np.zeros((self.naxis1, self.naxis2), dtype=np.float32)
        
        # At this point, the only negatives come from the convolution.
        # Remove them for the Poisson process
        self.array[self.array < 0] = 0
        
        out_array = np.zeros(self.array.shape, dtype=np.float32)
        
        if cmds["SIM_SPEED"] <= 3:
            for n in range(ndit):
                out_array += self._read_out_uptheramp(self.array, dit, tro,
                                                   max_byte)
                out_array += dark * dit

        elif cmds["SIM_SPEED"] > 3 and cmds["SIM_SPEED"] <= 7:
            for n in range(ndit):
                out_array += self._read_out_fast(self.array, dit)
                out_array += dark * dit
                
        elif cmds["SIM_SPEED"] > 7:
            out_array = self._read_out_superfast(self.array, dit, ndit)
            out_array += dark * dit * ndit
            
        #### TODO #########
        # add read out noise for every readout 
        out_array += self._read_noise_frame(cmds) * ndit

        if cmds["OBS_REMOVE_CONST_BG"].lower() == "yes":
            min_val = np.min(out_array)
            out_array -= min_val
        
        return out_array

    ## TODO: What to do if dit = min_dit (single read)?
    ## TODO: Make breaking up into memory chunks more flexible?
    def _read_out_uptheramp(self, image, dit, tro=1.3, max_byte=2**30):
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


        Keyword Arguments (**kwargs)
        ----------------------------


        Raises
        ------


        See Also
        --------


        Notes
        -----


        References
        ----------


        Examples
        --------
        """
        """Test readout onto a detector using cube model

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
        tpts =  (1 + np.arange(nro)) * tro

        img_byte = image.nbytes
        pix_byte = img_byte / (nx * ny)

        #max_byte =            ## TODO: arbitrary, function parameter?
        max_pix = max_byte / pix_byte

        cube_megabyte = img_byte * nro / 2**20
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
                    del(cube)
                    cube = np.zeros((nx, ny_cut, nro))
                except:
                    pass

            ## Fill the cube with Poisson realization, individual reads
            for i in range(nro):
                cube[:,:,i] = np.random.poisson(image[:,y1:y2] * tro)

            ## Build the ramp
            sumcube = cube.cumsum(axis=2)

            ## determine the slope using explicit formula calculated over cube
            Sx = tpts.sum()
            Sxx = (tpts * tpts).sum()
            Sy = np.sum(sumcube, axis=2)
            Sxy = np.sum(sumcube * tpts, axis=2)

            slope[:,y1:y2] = (nro * Sxy - Sx * Sy) / (nro * Sxx - Sx * Sx)

            ## Move to next slice
            y1 = y2

        # return values are [ph/pixel]
        return slope * dit

    def _read_out_fast(self, image, dit):
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


        Keyword Arguments (**kwargs)
        ----------------------------


        Raises
        ------


        See Also
        --------


        Notes
        -----


        References
        ----------


        Examples
        --------
        """

        return np.random.poisson(image * dit)

        
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


        Keyword Arguments (**kwargs)
        ----------------------------


        Raises
        ------


        See Also
        --------


        Notes
        -----


        References
        ----------


        Examples
        --------
        """

        # As noise increases with sqrt(t), we 
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
        
        if cmds["FPA_NOISE_PATH"] is not None:
            n = len(fits.info(cmds["FPA_NOISE_PATH"], False))
            layer = np.random.randint(n)
            tmp = fits.getdata(cmds["FPA_NOISE_PATH"], layer)
            return tmp[:self.naxis1, :self.naxis2]
        elif "gen" in cmds["FPA_NOISE_PATH"].lower():
            if cmds["HXRG_OUTPUT_PATH"] is not None:
                self._generate_hxrg_noise(cmds)
                tmp = fits.getdata(cmds["HXRG_OUTPUT_PATH"])
                return tmp[:self.naxis1, :self.naxis2]
            else:
                return self._generate_hxrg_noise(cmds)
        else:
            return np.zeros((self.naxis1, self.naxis2))
        
        
    def _generate_hxrg_noise(self, cmds):
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


        Keyword Arguments (**kwargs)
        ----------------------------


        Raises
        ------


        See Also
        --------


        Notes
        -----


        References
        ----------


        Examples
        --------
        """


        """
        Create a detector noise array using Bernard Rauscher's NGHxRG tool

        Optional Keywords:
        - HXRG_OUTPUT_PATH:

        """
        if len(kwargs) > 0 and self.verbose: print("updating ",kwargs)
        self.params.update(kwargs)

        # HXRG needs a pca file to run. Work out what a PCA file means!!
        ng_h4rg     = HXRGNoise(naxis1   = cmds["CHIP_NAXIS1"],
                                naxis2   = cmds["CHIP_NAXIS2"],
                                n_out    = cmds["HXRG_NUM_OUTPUTS"],
                                nroh     = cmds["HXRG_NUM_ROW_OH"],
                                pca0_file= cmds["HXRG_PCA0_FILENAME"],
                                verbose  = cmds["SIM_VERBOSE"])

        # Make a noise file
        noise = ng_h4rg.mknoise(o_file   = cmds["HXRG_OUTPUT_PATH"],
                                rd_noise = cmds["FPA_READOUT_MEDIAN"],
                                pedestal = cmds["HXRG_PEDESTAL"],
                                c_pink   = cmds["HXRG_CORR_PINK"],
                                u_pink   = cmds["HXRG_UNCORR_PINK"],
                                acn      = cmds["HXRG_ALT_COL_NOISE"])
        
        return noise


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
        
        
###############################################################################
#                       NGHXRG by Bernard Rauscher                            #
#             see the paper: http://arxiv.org/abs/1509.06264                  #
#           downloaded from: http://jwst.nasa.gov/publications.html           #
###############################################################################

# dependencies include: astropy, numpy, scipy [, datetime, warnings, os]
# import os
# import warnings
# from astropy.io import fits
# import numpy as np
# from scipy.ndimage.interpolation import zoom
# import datetime
# import matplotlib.pyplot as plt # Handy for debugging

#warnings.filterwarnings('ignore')

class HXRGNoise:
    """
    HXRGNoise is a class for making realistic Teledyne HxRG system
    noise. The noise model includes correlated, uncorrelated,
    stationary, and non-stationary components. The default parameters
    make noise that resembles Channel 1 of JWST NIRSpec. NIRSpec uses
    H2RG detectors. They are read out using four video outputs at
    1.e+5 pix/s/output.
    """

    # These class variables are common to all HxRG detectors
    nghxrg_version = 2.3 # Software version

    def __init__(self, naxis1=None, naxis2=None, naxis3=None, n_out=None,
                 dt=None, nroh=None, nfoh=None, pca0_file=None, verbose=False,
                 reverse_scan_direction=False,
                 reference_pixel_border_width=None):
        """
        Simulate Teledyne HxRG+SIDECAR ASIC system noise.

        Parameters:
            naxis1      - X-dimension of the FITS cube
            naxis2      - Y-dimension of the FITS cube
            naxis3      - Z-dimension of the FITS cube
                          (number of up-the-ramp samples)
            n_out       - Number of detector outputs
            nfoh        - New frame overhead in rows. This allows for a short
                          wait at the end of a frame before starting the next
                          one.
            nroh        - New row overhead in pixels. This allows for a short
                          wait at the end of a row before starting the next one.
            dt          - Pixel dwell time in seconds
            pca0_file   - Name of a FITS file that contains PCA-zero
            verbose     - Enable this to provide status reporting
            reference_pixel_border_width - Width of reference pixel border
                                           around image area
            reverse_scan_direction - Enable this to reverse the fast scanner
                                     readout directions. This
                                     capability was added to support
                                     Teledyne's programmable fast scan
                                     readout directions. The default
                                     setting =False corresponds to
                                     what HxRG detectors default to
                                     upon power up.
        """

        # ======================================================================
        #
        # DEFAULT CLOCKING PARAMETERS
        #
        # The following parameters define the default HxRG clocking pattern. The
        # parameters that define the default noise model are defined in the
        # mknoise() method.
        #
        # ======================================================================

        # Default clocking pattern is JWST NIRSpec
        self.naxis1    = 2048  if naxis1   is None else int(naxis1)
        self.naxis2    = 2048  if naxis2   is None else int(naxis2)
        self.naxis3    = 1     if naxis3   is None else int(naxis3)
        self.n_out     = 4     if n_out    is None else int(n_out)
        self.dt        = 1.e-5 if dt       is None else dt
        self.nroh      = 12    if nroh     is None else int(nroh)
        self.nfoh      = 1     if nfoh     is None else int(nfoh)
        self.reference_pixel_border_width = 4 \
                                            if reference_pixel_border_width is \
                                            None else reference_pixel_border_width

        # Initialize PCA-zero file and make sure that it exists and is a file
        self.pca0_file = os.getenv('NGHXRG_HOME')+'/nirspec_pca0.fits' if \
                         pca0_file is None else pca0_file
        if os.path.isfile(self.pca0_file) is False:
            print('There was an error finding pca0_file! Check to be')
            print('sure that the NGHXRG_HOME shell environment')
            print('variable is set correctly and that the')
            print('$NGHXRG_HOME/ directory contains the desired PCA0')
            print('file. The default is nirspec_pca0.fits.')
            os.sys.exit()


        # ======================================================================

        # Configure status reporting
        self.verbose = verbose

        # Configure readout direction
        self.reverse_scan_direction = reverse_scan_direction

        # Compute the number of pixels in the fast-scan direction per
        # output
        self.xsize = self.naxis1 // self.n_out

        # Compute the number of time steps per integration, per
        # output
        self.nstep = (self.xsize+self.nroh) * (self.naxis2+self.nfoh)\
                     * self.naxis3

        # For adding in ACN, it is handy to have masks of the even
        # and odd pixels on one output neglecting any gaps
        self.m_even = np.zeros((self.naxis3,self.naxis2,self.xsize))
        self.m_odd = np.zeros_like(self.m_even)
        for x in np.arange(0,self.xsize,2):
            self.m_even[:,:self.naxis2,x] = 1
            self.m_odd[:,:self.naxis2,x+1] = 1
        self.m_even = np.reshape(self.m_even, np.size(self.m_even))
        self.m_odd = np.reshape(self.m_odd, np.size(self.m_odd))

        # Also for adding in ACN, we need a mask that point to just
        # the real pixels in ordered vectors of just the even or odd
        # pixels
        self.m_short = np.zeros((self.naxis3, self.naxis2+self.nfoh, \
                                      (self.xsize+self.nroh)//2))
        self.m_short[:,:self.naxis2,:self.xsize//2] = 1
        self.m_short = np.reshape(self.m_short, np.size(self.m_short))

        # Define frequency arrays
        self.f1 = np.fft.rfftfreq(self.nstep) # Frequencies for nstep elements
        self.f2 = np.fft.rfftfreq(2*self.nstep) # ... for 2*nstep elements

        # Define pinkening filters. F1 and p_filter1 are used to
        # generate ACN. F2 and p_filter2 are used to generate 1/f noise.
        self.alpha = -1 # Hard code for 1/f noise until proven otherwise
        self.p_filter1 = np.sqrt(self.f1**self.alpha)
        self.p_filter2 = np.sqrt(self.f2**self.alpha)
        self.p_filter1[0] = 0.
        self.p_filter2[0] = 0.


        # Initialize pca0. This includes scaling to the correct size,
        # zero offsetting, and renormalization. We use robust statistics
        # because pca0 is real data
        hdu = fits.open(self.pca0_file)
        naxis1 = hdu[0].header['naxis1']
        naxis2 = hdu[0].header['naxis2']
        if (naxis1 != self.naxis1 or naxis2 != self.naxis2):
            zoom_factor = self.naxis1 / naxis1
            self.pca0 = zoom(hdu[0].data, zoom_factor, order=1, mode='wrap')
        else:
            self.pca0 = hdu[0].data
        self.pca0 -= np.median(self.pca0) # Zero offset
        self.pca0 /= (1.4826*mad(self.pca0)) # Renormalize


    def message(self, message_text):
        """
        Used for status reporting
        """
        if self.verbose is True:
            print('NG: ' + message_text + ' at DATETIME = ', \
                  datetime.datetime.now().time())

    def white_noise(self, nstep=None):
        """
        Generate white noise for an HxRG including all time steps
        (actual pixels and overheads).

        Parameters:
            nstep - Length of vector returned
        """
        return(np.random.standard_normal(nstep))

    def pink_noise(self, mode):
        """
        Generate a vector of non-periodic pink noise.

        Parameters:
            mode - Selected from {'pink', 'acn'}
        """

        # Configure depending on mode setting
        if mode is 'pink':
            nstep = 2*self.nstep
            f = self.f2
            p_filter = self.p_filter2
        else:
            nstep = self.nstep
            f = self.f1
            p_filter = self.p_filter1

        # Generate seed noise
        mynoise = self.white_noise(nstep)

        # Save the mean and standard deviation of the first
        # half. These are restored later. We do not subtract the mean
        # here. This happens when we multiply the FFT by the pinkening
        # filter which has no power at f=0.
        the_mean = np.mean(mynoise[:nstep//2])
        the_std = np.std(mynoise[:nstep//2])

        # Apply the pinkening filter.
        thefft = np.fft.rfft(mynoise)
        thefft = np.multiply(thefft, p_filter)
        result = np.fft.irfft(thefft)
        result = result[:nstep//2] # Keep 1st half

        # Restore the mean and standard deviation
        result *= the_std / np.std(result)
        result = result - np.mean(result) + the_mean

        # Done
        return(result)



    def mknoise(self, o_file, rd_noise=None, pedestal=None, c_pink=None,
                u_pink=None, acn=None, pca0_amp=None,
                reference_pixel_noise_ratio=None, ktc_noise=None,
                bias_offset=None, bias_amp=None):
        """
        Generate a FITS cube containing only noise.

        Parameters:
            o_file   - Output filename
            pedestal - Magnitude of pedestal drift in electrons
            rd_noise - Standard deviation of read noise in electrons
            c_pink   - Standard deviation of correlated pink noise in electrons
            u_pink   - Standard deviation of uncorrelated pink noise in
                       electrons
            acn      - Standard deviation of alterating column noise in
                       electrons
            pca0     - Standard deviation of pca0 in electrons
            reference_pixel_noise_ratio - Ratio of the standard deviation of
                                          the reference pixels to the regular
                                          pixels. Reference pixels are usually
                                          a little lower noise.
            ktc_noise   - kTC noise in electrons. Set this equal to
                          sqrt(k*T*C_pixel)/q_e, where k is Boltzmann's
                          constant, T is detector temperature, and C_pixel is
                          pixel capacitance. For an H2RG, the pixel capacitance
                          is typically about 40 fF.
            bias_offset - On average, integrations stare here in electrons. Set
                          this so that all pixels are in range.
            bias_amp    - A multiplicative factor that we multiply PCA-zero by
                          to simulate a bias pattern. This is completely
                          independent from adding in "picture frame" noise.

        Note1:
        Because of the noise correlations, there is no simple way to
        predict the noise of the simulated images. However, to a
        crude first approximation, these components add in
        quadrature.

        Note2:
        The units in the above are mostly "electrons". This follows convention
        in the astronomical community. From a physics perspective, holes are
        actually the physical entity that is collected in Teledyne's p-on-n
        (p-type implants in n-type bulk) HgCdTe architecture.
        """

        self.message('Starting mknoise()')

        # ======================================================================
        #
        # DEFAULT NOISE PARAMETERS
        #
        # These defaults create noise similar to that seen in the JWST NIRSpec.
        #
        # ======================================================================

        self.rd_noise  = 5.2      if rd_noise     is None else rd_noise
        self.pedestal  = 4        if pedestal     is None else pedestal
        self.c_pink    = 3        if c_pink       is None else c_pink
        self.u_pink    = 1        if u_pink       is None else u_pink
        self.acn       = .5       if acn          is None else acn
        self.pca0_amp  = .2       if pca0_amp     is None else pca0_amp

        # Change this only if you know that your detector is different from a
        # typical H2RG.
        self.reference_pixel_noise_ratio = 0.8 if \
            reference_pixel_noise_ratio is None else reference_pixel_noise_ratio

        # These are used only when generating cubes. They are
        # completely removed when the data are calibrated to
        # correlated double sampling or slope images. We include
        # them in here to make more realistic looking raw cubes.
        self.ktc_noise   = 29.   if ktc_noise   is None else ktc_noise
        self.bias_offset = 5000. if bias_offset is None else bias_offset
        self.bias_amp    = 500.  if bias_amp    is None else bias_amp

        # ======================================================================

        # Initialize the result cube. For up-the-ramp integrations,
        # we also add a bias pattern. Otherwise, we assume
        # that the aim was to simulate a two dimensional correlated
        # double sampling image or slope image.
        self.message('Initializing results cube')
        result = np.zeros((self.naxis3, self.naxis2, self.naxis1), \
                          dtype=np.float32)
        if self.naxis3 > 1:
            # Inject a bias pattern and kTC noise. If there are no reference pixels,
            # we know that we are dealing with a subarray. In this case, we do not
            # inject any bias pattern for now.
            if self.reference_pixel_border_width > 0:
                bias_pattern = self.pca0*self.bias_amp + self.bias_offset
            else:
                bias_pattern = self.bias_offset

            # Add in some kTC noise. Since this should always come out
            # in calibration, we do not attempt to model it in detail.
            bias_pattern += \
                         self.ktc_noise * \
                         np.random.standard_normal((self.naxis2, self.naxis1))

            # Ensure that there are no negative pixel values. Data cubes
            # are converted to unsigned integer before writing.
            bias_pattern = np.where(bias_pattern < 0, 0, bias_pattern)

            # Add in the bias pattern
            for z in np.arange(self.naxis3):
                result[z,:,:] += bias_pattern


        # Make white read noise. This is the same for all pixels.
        self.message('Generating rd_noise')
        w = self.reference_pixel_border_width # Easier to work with
        r = self.reference_pixel_noise_ratio  # Easier to work with
        for z in np.arange(self.naxis3):
            here = np.zeros((self.naxis2, self.naxis1))
            if w > 0: # Ref. pixel border exists
                # Add both reference and regular pixels
                here[:w,:] = r * self.rd_noise * \
                             np.random.standard_normal((w,self.naxis1))
                here[-w:,:] = r * self.rd_noise * \
                              np.random.standard_normal((w,self.naxis1))
                here[:,:w] = r * self.rd_noise * \
                             np.random.standard_normal((self.naxis1,w))
                here[:,-w:] = r * self.rd_noise * \
                              np.random.standard_normal((self.naxis1,w))
                # Make noisy regular pixels
                here[w:-w,w:-w] = self.rd_noise * \
                                  np.random.standard_normal( \
                                  (self.naxis2-2*w,self.naxis1-2*w))
            else: # Ref. pixel border does not exist
                # Add only regular pixels
                here = self.rd_noise * np.random.standard_normal((self.naxis2,\
                                                                  self.naxis1))
            # Add the noise in to the result
            result[z,:,:] += here


        # Add correlated pink noise.
        self.message('Adding c_pink noise')
        tt = self.c_pink * self.pink_noise('pink') # tt is a temp. variable
        tt = np.reshape(tt, (self.naxis3, self.naxis2+self.nfoh, \
                             self.xsize+self.nroh))[:,:self.naxis2,:self.xsize]
        for op in np.arange(self.n_out):
            x0 = op * self.xsize
            x1 = x0 + self.xsize
            if self.reverse_scan_direction is False:
                # Teledyne's default fast-scan directions
                if np.mod(op,2)==0:
                    result[:,:,x0:x1] += tt
                else:
                    result[:,:,x0:x1] += tt[:,:,::-1]
            else:
                # Reverse the fast-scan directions.
                if np.mod(op,2)==1:
                    result[:,:,x0:x1] += tt
                else:
                    result[:,:,x0:x1] += tt[:,:,::-1]



        # Add uncorrelated pink noise. Because this pink noise is stationary and
        # different for each output, we don't need to flip it.
        self.message('Adding u_pink noise')
        for op in np.arange(self.n_out):
            x0 = op * self.xsize
            x1 = x0 + self.xsize
            tt = self.u_pink * self.pink_noise('pink')
            tt = np.reshape(tt, (self.naxis3, self.naxis2+self.nfoh, \
                             self.xsize+self.nroh))[:,:self.naxis2,:self.xsize]
            result[:,:,x0:x1] += tt

        # Add ACN
        self.message('Adding acn noise')
        for op in np.arange(self.n_out):

            # Generate new pink noise for each even and odd vector.
            # We give these the abstract names 'a' and 'b' so that we
            # can use a previously worked out formula to turn them
            # back into an image section.
            a = self.acn * self.pink_noise('acn')
            b = self.acn * self.pink_noise('acn')

            # Pick out just the real pixels (i.e. ignore the gaps)
            a = a[np.where(self.m_short == 1)]
            b = b[np.where(self.m_short == 1)]

            # Reformat into an image section. This uses the formula
            # mentioned above.
            acn_cube = np.reshape(np.transpose(np.vstack((a,b))),
                                  (self.naxis3,self.naxis2,self.xsize))

            # Add in the ACN. Because pink noise is stationary, we can
            # ignore the readout directions. There is no need to flip
            # acn_cube before adding it in.
            x0 = op * self.xsize
            x1 = x0 + self.xsize
            result[:,:,x0:x1] += acn_cube


        # Add PCA-zero. The PCA-zero template is modulated by 1/f.
        if self.pca0_amp > 0:
            self.message('Adding PCA-zero "picture frame" noise')
            gamma = self.pink_noise(mode='pink')
            zoom_factor = self.naxis2 * self.naxis3 / np.size(gamma)
            gamma = zoom(gamma, zoom_factor, order=1, mode='mirror')
            gamma = np.reshape(gamma, (self.naxis3,self.naxis2))
            for z in np.arange(self.naxis3):
                for y in np.arange(self.naxis2):
                    result[z,y,:] += self.pca0_amp*self.pca0[y,:]*gamma[z,y]


        # If the data cube has only 1 frame, reformat into a 2-dimensional
        # image.
        if self.naxis3 == 1:
            self.message('Reformatting cube into image')
            result = result[0,:,:]

        # If the data cube has more than one frame, convert to unsigned
        # integer
        if self.naxis3 > 1:
            self.message('Converting to 16-bit unsigned integer')
            result = result.astype('uint16')

        # Write the result to a FITS file
        self.message('Writing FITS file')
        hdu = fits.PrimaryHDU(result)
        hdu.header.append()
        hdu.header.append(('RD_NOISE', self.rd_noise, 'Read noise'))
        hdu.header.append(('PEDESTAL', self.pedestal, 'Pedestal drifts'))
        hdu.header.append(('C_PINK', self.c_pink, 'Correlated pink'))
        hdu.header.append(('U_PINK', self.u_pink, 'Uncorrelated pink'))
        hdu.header.append(('ACN', self.acn, 'Alternating column noise'))
        hdu.header.append(('PCA0', self.pca0_amp, \
                           'PCA zero, AKA picture frame'))
        #hdu.header['HISTORY'] = 'Created_by_NGHXRG_version_' \
        #                        + str(self.nghxrg_version)

        self.message('Exiting mknoise()')

        if o_file is not None:
            hdu.writeto(o_file, clobber='True')
        return result
        
class bloedsinn:
    pass