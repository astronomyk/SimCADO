########################################################################
#                          SIMCADO v0.0                                #
#             Created 21.11.2015 by Kieran Leschinski                  #
########################################################################

# This is a shell document to collect the Objects and Methods needed
# for the the alpha release of SimCADO

# Coding Style (so that Kieran sticks to it):
# my_variable_name
# myFunctionName() 
# MyClassName


########################################################################
#    UserCommands 
########################################################################

#class UserCommands:
"""
Read in the configuration files
The UserCommands can be held in a configuration file or generates by 
the user interactively. The commands then generate two objects
- AtmosphericModel
- LightObject
"""

    #def __init__(self, filename):
    """
    """

########################################################################
#    AtmosphericModel
########################################################################

#class AtmosphericModel:
"""
Info about the Atmosphere
"""

    #def __init__(from_file=True):
    """
    Load in a SkyCalc FITS file or 
    """


########################################################################
#    LightObject
########################################################################

#class LightObject:

"""
Info about the Telescope and the Instrument
"""

    #def __init__():
    """
    """
	
    #def add_atmosphere()
	
    #def apply_psf()
	
    #def apply_throughput()
	
    #def collapse_cube()
	
    #def apply_plane_effects()
	
    #def make_detector_image()
	
    #def add_background_light()
	


########################################################################
#    OpticalTrain
########################################################################

#class OpticalTrain:

"""
Info about the Telescope and the Instrument

The OpticalTrain class creates several subclasses
- PsfCube
- ThroughputCurve
- PlaneEffect
- Detector

Any of these objects can be created by the user and replaced in the main
instance of the OpticalTrain 

"""

    #def __init__():
    """
    """

    #def make_psf_cube()
	
    #def make_throughput_curve()
	
    #def make_plane_effect()
	
    #def make_detector
	
########################################################################
#    PSFCube
########################################################################

#class PSFCube():
"""
Info about the Telescope and the Instrument
"""

    #def __init__(self, optical_train = None):
    """
    """
	
########################################################################
#    ThroughputCurve
########################################################################

#class ThroughputCurve:
"""
Generate a 
"""

    #def __init__(self, optical_train = None, kwargs**):
    """
    Import all the spectral curves that are listed in the optical_train
    or those in the kwargs (filenames, labels)
    
    List of self variables
    self.filenames
    self.labels
    self.lam
    self.val
    self.dlam
    """

    #def collapse(self, labels)
    """
    Combine all the curves specified in "labels". That way we can change the 
    throughput for each source of photons:
    - science object - atmo + M1->M7 + window + filter + QE
    - M1 greybody emission - M2->M7 + window + filter + QE
    """
    
	
########################################################################
#    PlaneEffect
########################################################################

#class PlaneEffect:
"""
Generate whatever we need to apply the plane effects
- rotation
- tracking errors
- wind shake
- distortion
"""

    #def __init__(self, optical_train = None):
    """
    if optical_train is not None:
    - import the distortion map, if there is one
    - generate a PSF based on the wind shake
    - generate a PSF for the tracking errors
    - generate the derotator error (whatever this entails)
    else:
    - create a list of empty variables
    - wind_shake_psf = None
    - tracking_error_psf = None
    - derotator_effectness = None
    
    List of self variables
    self.distortion_map [new type?]
    self.wind_shake_psf [PSF type - 2D Gaussian Elliptical]
    self.tracking_error_psf [PSF type - 1D Gaussian]
    
    """
    
    
########################################################################
#    Detector
########################################################################

#class Detector:
"""
- check to see if there is a file name for the distortion map in the 
    optical_train
- generate the detector characteristics
    - everything from NG_HxRG?
"""

    #def __init__(self, optical_train=None, kwargs**):
    """
    Leave this for Oliver
    # use NG_HxRG to generate the various arrays needed?
    
    List of self variables (???)
    # all arrays
    self.read_noise [2D array]
    self.dark_noise [2D array]
    self.pixel_sensitivities [2D array]
    ...
    """


########################################################################
#    SpectralCurve and PSFLayer
########################################################################

#class SpectralCurve:
"""
Info about the a spectral curve element
"""

    #def __init__(self, from_file = True, lam = None, lam_unit = None, 
    #             val = None, val_unit = None, label = None):
    """
    If from_file == True: load in the file -- test for ascii or fits
    Else: use the values given by the user
    
    List of self variables
    self.lam [1D array]
    self.val [1D array]
    self.dlam [float]
    self.lam_orig [list of 1D array]
    self.val_orig [list of 1D array]
    self.label [string]
    """
    
    #def add(self, spectral_curve, label = None)
    """
    Add another spectral curve to the list
    Is it a function or a method?
    """

    #def collapse(self, dlam = None)
    """
    Take all the curves in self.lam_orig and self.val_orig and combine 
    them into a single curve at a resolution specified by dlam
    Put the output into self.lam, self.val, self.dlam
    """
    
    #def rescale(self, dlam = 1*u.nm)
    """
    Rescale ONLY self.lam, self.val to the dlam bin width
    Update self.dlam
    """
    
    
#class PSFLayer    
"""
The PSFLayer is just a simple array with dimensions (2^n, 2^n) [where
2^n is the next biggest power of 2 for the needed size of the PSF array
"""    
    #def __init__(self, type, kwargs**)
    """
    Types of PSFs (and their needed keyword arguments):
    - Gaussian (elliptical)
        kwargs - (x,y,fwhm,ellipticity=0,angle=0)
    - Airy
        kwargs - (x,y,fwhm)
    - Moffat
        kwargs - (???)
    - Dot
        kwargs - (x,y)
    - Strip
        kwargs - (x,y,angle,vector)
    - 1D Gaussian at an angle
        kwargs - (x,y,angle,fwhm)
    - from FITS file
        kwargs - (filename)
    - array from user
    
    List of self variables
    self.psf
    
    """





