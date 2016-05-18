"""
This module contains classes which control how a simulation is run


Summary
-------
UserCommands is essentially a dictionary that holds all the variables that
the user may wish to change. It also has some set variables like `pix_res`
that can be accessed directly, instead of from the dictionary.

UserCommands is imported directly into the simcado package and is accessable 
from the main package - `simcado.UserCommands`

If UserCommands is called without any arguments, the default values for MICADO 
are used. 


Routines
--------
`UserCommands(filename, default=<path_to_default>)`


See Also
--------
Class that require a `UserCommands` object directly include:
- `Detector`
- `OpticalTrain`


Notes
-----

References
----------

Examples
--------
By default `UserCommands` contains the parameters needed to generate the MICADO
optical train:
```
>>> import simcado
>>> my_cmds = simcado.UserCommands()
>>> my_cmds["SCOPE_M1_DIAMETER_OUT"]
37.3
```

To list the keywords that are available:
```
>>> my_cmds.keys()
...
```

The UserCommands object also contains smaller dictionaries for each category of
keywords - e.g. for the keywords for the instrument:
```
>>> my_cmds.inst
...
```

"""


import os, warnings
import numpy as np
try:
    import simcado.spectral as sc
    import simcado.utils as utils
except:
    import spectral as sc
    import utils as utils

    
__all__ = ["UserCommands"]    
    

class UserCommands(object):
    """
    An extended dictionary with the parameters needed for running a simulation

    
    Summary
    -------
    A `UserCommands` object contains a dictionary which holds all the keywords 
    from the `default.config` file. It also has attributes which represent the
    frequently used variables, i.e. `pix_res`, `lam_bin_edges`, `exptime`, etc
    
    `<UserCommands>.cmds` is a dictionary that holds all the variables
    the user may wish to change. It also has some set variables like 
    `<UserCommands>.pix_res` that can be accessed directly, instead of from the 
    dictionary.

    `UserCommands` is imported directly into the simcado package and is 
    accessable from the main package - `simcado.UserCommands`

    If UserCommands is called without any arguments, the default values for 
    MICADO and the E-ELT are used. 

    
    Parameters
    ----------
    filename : str
        path to the user's .config file
    default : str, optional
        path to the default.config file.
    
    
    Attributes
    ----------
    Internal dictionaries
    cmds : dict
        the dictionary which holds all the keyword-value pairs needed for 
        running a simualtion
    obs : dict
        parameters about the observation
    sim : dict
        parameters about the simualtion
    atmo : dict
        parameters about the atmosphere
    scope : dict
        parameters about the telescope
    inst : dict
        parameters about the instrument
    fpa : dict
        parameters about the detector array (FPA - Focal Plane Array)
    hxrg : dict
        parameters about the chip noise (HxRG - HAWAII 4RG chip series)
    
    Attributes pertaining to the purely spectral data sets (e.g. transmission
    curves, stellar spectra, etc)
    lam : np.ndarray
        a vector containing the centres of the wavelength bins used when 
        resampling the spectra or transmission curves
    lam_res : float
        [um] the resolution of the `lam`
        
    Attributes pertaining to the binning in spectral space for which different
    PSFs need to be used
    lam_psf_res : float
        [um] the spectal "distance" between layers - i.e. width of the bins
    lam_bin_edges : array-like
        [um] the edge of the spectral bin for each layer
    lam_bin_centers : array-like
        [um] the centres of the spectral bin for each layer
    
    Attributes pertaining to the binning in the spatial plane
    pix_res : float
        [arcsec] the internal (oversampled) spatial resolution of the simulation
    fpa_res : float
        [arcsec] the field of view of the individual pixels on the detector
    
    General attributes
    verbose : bool
        Flag for printing intermediate results to the screen (default=True)
    exptime : float
        [s] exposure time of a single DIT
    diameter : float
        [m] outer diamter of the primary aperture (i.e. M1)
    area : float
        [m^2] effective area of the primary aperture (i.e. M1)
    
    
    Methods
    -------
    update(new_dict)
        updates the current `UserCommands` object from another dict-like object
    keys()
        returns the keys in the `UserCommands.cmds` dictionary
    values()
        returns the values in the `UserCommands.cmds` dictionary
    
    Raises
    ------

    See Also
    --------
    Detector, OpticalTrain

    Notes
    -----

    References
    ----------

    Examples
    --------
    By default `UserCommands` contains the parameters needed to generate the 
    MICADO optical train:
    ```
    >>> import simcado
    >>> my_cmds = simcado.UserCommands()
    >>> my_cmds["SCOPE_M1_DIAMETER_OUT"]
    37.3
    ```
    
    `UserCommands` supports indexing like a dictionary object.
    ```
    >>> my_cmds["SCOPE_M1_DIAMETER_OUT"] = 8.2
    >>> my_cmds["SCOPE_M1_DIAMETER_OUT"]
    8.2
    ```

    To list the keywords that are available:
    ```
    >>> my_cmds.keys()
    ...
    ```

    The `UserCommands` object also contains smaller dictionaries for each category
    of keywords - e.g. for the keywords describing the instrument:
    ```
    >>> my_cmds.inst
    ...
    ```
    """
    
    
    def __init__(self, filename=None,
                 default="../data/default.config"):

        """
        Create an extended dictionary of simulation parameters
        """
         
        self.cmds = utils.read_config(default)

        if user_filename is not None:
            self.cmds.update(utils.read_config(filename))

        self.cmds["CONFIG_USER"]   = user_filename
        self.cmds["CONFIG_DEFAULT"] = default_filename

        # check the output path directory and file name
        if self.cmds["OBS_OUTPUT_DIR"] is None:
            self.cmds["OBS_OUTPUT_DIR"] = "./"

        if self.cmds["OBS_OUTPUT_NAME"] is None:
            self.cmds["OBS_OUTPUT_NAME"] = "output.fits"

        if self.cmds["SIM_PSF_OVERSAMPLE"] == "yes":
            self.cmds["PSF_MODE"] = "oversample"
        else:
            self.cmds["PSF_MODE"] = "linear_interp"

        # Check for a filter curve file or a standard broadband name
        if self.cmds["INST_FILTER_TC"] in ["I", "z", "Y", "J", "H", "Ks", "K"]:
            self.cmds["INST_FILTER_TC"] = "../data/TC_filter_" + \
                                            self.cmds["INST_FILTER_TC"] + ".dat"

        self._update_attributes()

        if self.verbose:
            print("Read in parameters from \n"+"\n".join(fnames))

            
    def update(self, new_dict):
        """
        Update nultiple entries of a `UserCommands` dictionary


        Summary
        -------
        `update(new_dict)` takes either a normal python `dict` object or a 
        `UserCommands` object. Only keywords that match those in the 
        `UserCommands` object will be updated. No warning is given for 
        misspelled keywords.
        
        To update single items in the dictionary, it is recommended to simply
        call the key and update the value - i.e `<UserCommands>[key] = value`. 
        
        Parameters
        ----------
        new_dict : dict, `UserCommands`

       
        Raises
        ------

        See Also
        --------
        UserCommands
        
        Notes
        -----

        Examples
        --------
        ```
        >>> import simcado
        >>> my_cmds = simcado.UserCommands()
        >>> new_cmds = {"OBS_EXPTIME" : 60 , "OBS_NDIT" : 10}
        >>> my_cmds.update(new_cmds)
        ```
        """

        if isinstance(new_dict, commands):
            self.cmds.update(new_dict.cmds)
        elif isinstance(new_dict, dict):
            self.cmds.update(new_dict)
        else:
            raise ValueError("Cannot update with type: "+type(new_dict))
        self._update_attributes()

        
    def keys(self):
        """
        Return the keys in the `UserCommands.cmds` dictionary
        """
        return self.cmds.keys()

    def values(self):
        """
        Return the values in the `UserCommands.cmds` dictionary
        """
        return self.cmds.values()
            
            
            

    def _update_attributes(self):
        """
        Update the UserCommand convenience attributes
        """

        self.fpa_res    = self.cmds["SIM_DETECTOR_PIX_SCALE"]
        self.pix_res    = self.fpa_res / self.cmds["SIM_OVERSAMPLING"]

        # if SIM_USE_FILTER_LAM is true, then use the filter curve to set the
        # wavelength boundaries where the filter is < SIM_FILTER_THRESHOLD

        if self.cmds["SIM_USE_FILTER_LAM"].lower() == "yes":
            tc_filt = sc.TransmissionCurve(filename=self.cmds['INST_FILTER_TC'])
            mask = np.where(tc_filt.val > self.cmds["SIM_FILTER_THRESHOLD"])[0]
            i0 = np.max((mask[0] - 1, 0))
            i1 = np.min((mask[-1] + 1 , len(tc_filt.lam) - 1))
            lam_min, lam_max   = tc_filt.lam[i0], tc_filt.lam[i1]
        else:
            lam_min, lam_max = self.cmds["SIM_LAM_MIN"], self.cmds["SIM_LAM_MAX"]

        self.lam_res = self.cmds["SIM_LAM_TC_BIN_WIDTH"]
        self.lam = np.arange(lam_min, lam_max + 1E-7, self.lam_res)

        #self.lam_psf_res = self.cmds["SIM_LAM_PSF_BIN_WIDTH"]
        #self.lam_bin_edges = np.arange(lam_min,
        #                               lam_max + self.lam_psf_res + 1E-7,
        #                               self.lam_psf_res)
        # make lam_bin_edges according to how great the ADC offsets are
        self.lam_bin_edges = self._get_lam_bin_edges(lam_min, lam_max)
        self.lam_bin_centers = 0.5 * (self.lam_bin_edges[1:] + \
                                      self.lam_bin_edges[:-1])

        self.exptime    = self.cmds["OBS_EXPTIME"]
        self.diameter   = self.cmds["SCOPE_M1_DIAMETER_OUT"]
        self.area       = np.pi / 4 * (self.diameter**2 - \
                                        self.cmds["SCOPE_M1_DIAMETER_IN"]**2)

        self.cmds["SIM_N_MIRRORS"] = self.cmds["SCOPE_NUM_MIRRORS"] + \
                                     self.cmds["INST_NUM_MIRRORS"] + \
                                     self.cmds["INST_NUM_EXT_MIRRORS"]

        for key,value in zip(self.cmds.keys(), self.cmds.values()):
            if type(value) == str and value.lower() == "none":
                self.cmds[key] = None

        self.verbose = True   if self.cmds["SIM_VERBOSE"] == "yes"   else False

        self._split_categories()


    def _get_lam_bin_edges(self, lam_min, lam_max):
        """
        Generates an array with the bin edges of the layers in spectral space
        
        
        Parameters
        ----------
        lam_min, lam_max : float
            [um] the minimum and maximum wavelengths of the filter range
        
        
        Notes
        -------
        Atmospheric diffraction causes blurring in an image. To model this 
        effect the spectra from a `Source` object are cut into bins based on 
        how far the photons at different wavelength are diffracted from the 
        image center. The threshold for defining a new layer based on the how
        far a certain bin will move is given by `SIM_ADC_SHIFT_THRESHOLD`. The
        default value is 1 pixel.
        """
        
        effectiveness = self.cmds["INST_ADC_PERFORMANCE"] / 100.
        if effectiveness == 1.:
            lam_bin_edges = np.array([lam_min, lam_max])
            return lam_bin_edges
            
        shift_threshold = self.cmds["SIM_ADC_SHIFT_THRESHOLD"]

        ## get the angle shift for each slice
        lam = np.arange(lam_min, lam_max + 1E-7, 0.001)
        angle_shift = utils.atmospheric_refraction( lam,
                                                    self.cmds["OBS_ZENITH_DIST"],
                                                    self.cmds["ATMO_TEMPERATURE"],
                                                    self.cmds["ATMO_REL_HUMIDITY"],
                                                    self.cmds["ATMO_PRESSURE"],
                                                    self.cmds["SCOPE_LATITUDE"],
                                                    self.cmds["SCOPE_ALTITUDE"])

        ## convert angle shift into number of pixels
        ## pixel shifts are defined with respect to last slice
        rel_shift = (angle_shift - angle_shift[-1]) / self.pix_res
        if np.max(np.abs(rel_shift)) > 1000:
            raise ValueError("Pixel shifts too great (>1000), check units")

        ## Rotate by the paralytic angle
        rel_shift *= (1. - effectiveness)
        int_shift = np.array(rel_shift / shift_threshold, dtype=np.int)
        idx = [np.where(int_shift == i)[0][0]
               for i in np.unique(int_shift)[::-1]]
        lam_bin_edges = np.array(lam[idx + [len(lam)-1]])

        return lam_bin_edges


    def _split_categories(self):
        """
        Generate smaller category-specific dictionaries
        """
        self.obs    = {i:self.cmds[i] for i in self.cmds.keys() if "OBS" in i}
        self.sim    = {i:self.cmds[i] for i in self.cmds.keys() if "SIM" in i}
        self.atmo   = {i:self.cmds[i] for i in self.cmds.keys() if "ATMO" in i}
        self.scope  = {i:self.cmds[i] for i in self.cmds.keys() if "SCOPE" in i}
        self.inst   = {i:self.cmds[i] for i in self.cmds.keys() if "INST" in i}
        self.fpa    = {i:self.cmds[i] for i in self.cmds.keys() if "FPA" in i}
        self.hxrg   = {i:self.cmds[i] for i in self.cmds.keys() if "HXRG" in i}

   
    def __str__(self):
        if self.cmds["CONFIG_USER"] is not None:
            return "A dictionary of commands compiled from " + \
                                                        self.cmds["CONFIG_USER"]
        else:
            return "A dictionary of default commands"

    def __iter__(self):
        return self.cmds.__iter__()

    def __getitem__(self, key):
        return self.cmds[key]

    def __setitem__(self, key, val):
        self.cmds[key] = val
        self._update_attributes()


    ### Add to the update that all the cmds.variable are updated when
    ### the dicts are updated

    
class __bloedsinn():
    pass