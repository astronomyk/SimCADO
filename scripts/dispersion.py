#spectroscopy.py
"""
###############################
Projecting a slice onto the FPA
###############################


Each monochromatic slit image is projected to different coordinate(s) on the
focal plane.
    - class Slit
        - x_cen, y_cen, length, width [arcsec], angle_from_horizontal [degree]
    - SPEC_SLIT_ANGLE
    - SPEC_REF_WAVELENGTH

The position of a slice on the FPA is given by the trace_maps.
Read in all the trace_maps at the beginning and make a Trace class for each one

    - class Trace
        - get_focal_plane_positions(lambda)
              return [x],[y],[i]

    - SPEC_TRACE_LIST
        - trace_map_list.tbl
            - order, lam_min, lam_max, throughput, filename

Before the slit there is no ADC, thus the atmosphere stretches the PSF depending
on wavelength.
Calculate the offset of the slit w.r.t the FOV centre

    - OBS_ZENITH_DIST
    - ATMO_TEMPERATURE
    - ATMO_PRESSURE
    - ATMO_REL_HUMIDITY
    - .get_atmo_offset(lam, ref_lam, airmass, temp, pressure, rel_hum,
                                                            adc_effectivity=0)
        - returns x,y

As the orientation of the slit is fixed w.r.t. the FPA, the angle of the slit
results in an effective rotation of the source positions and the atmospheric
shifts.

    - src2 = copy(src)
    - src2.rotate(SPEC_SLIT_ANGLE)

Create a Chip which has the dimensions of the slit, but with the rotated x,y
offsets due to the slit angle and atmospheric dispersion.
The source should be imaged onto this Chip, but ONLY the wavelength range of the
current monochromatic slice should be used.

    - SIM_DETECTOR_PIX_SCALE
    - slit_chip = Chip(x_cen+dx_atmo, y_cen+dy_atmo, length, width, pix_res)
    - im = src2.image_in_range(lam, lam+dlam, slit_chip, psf)
    - <Slit>.apply_image(im)

Go through the list of Traces and get the position(s) on the focal plane where
the Trace needs to be splatted.
For each position, got through the list of Chips in the Detector to find onto
which detector they need to be spatted - may be up to 4 detectors, think corners

    - <Slit>.get_chips(<Detector>)
          return chip_list

Splat the slit image onto the require chips

    - <Slit>.splat_on_chips(chip_list)



##################################
Running a spectrograph observation
##################################

Each observation consists of creating a slit image of the <Source> object given
the atmospheric component (plus any other shifts?) for each wavelength bin
First steps
    - read in the trace_maps
    - find out the wavelengths that need to be projected. At the very least (but
        also the easiest) this is one image per pixel along a trace_map.
    - Go through the trace_maps and work out a list of wavelengths which correspond
        to a single pixel shift, at least in the spectral axis
        Sideways shift will need to be included (like astrometry, with a flag)

        <Disperser>.get_trace_map_list(SPEC_TRACE_LIST)
        <Disperser>.get_required_wavelengths(trace_map_list)

Make an <OpticalTrain>
Accept a <Source>
Apply the system transmission curve to the source
Add the background emission, essentially making a full flux spectrum

Loop over each of the wavelengths that need to be projected

Project slit image onto slice

Read out the chips into a fits file


####################
Extract the spectrum
####################

SimCADO should only make the detector images. However a script to extract the
spectra would be useful. Hence add an extra function.

The function should take find the wavelength for each of the integer pixel
values along the spectral axis of each trace_map. Then find the x,y positions of
the centre and extract and sum the values at position+/-width//2
All the trace_maps should then resampled and added together

    - extract_spectrum(hdu, trace_maps, position, width)
          return EmissionCurve or (lam, val)


"""

import os

import glob
import warnings
import logging
from copy import deepcopy

import numpy as np

from astropy.io import fits
from astropy.io import ascii as ioascii
from astropy.io.fits import BinTableHDU
from astropy.table import Table, Column
import astropy.units as u

from .optics import OpticalTrain
from .detector import Chip
from .commands import UserCommands
from .utils import atmospheric_refraction
from .utils import __pkg_dir__, __data_dir__

__all__ = []



class Trace(object):
    """
    Contains the parameters to map a slit image to a position on a detector chip


    Parameters
    ----------
    x : astropy.BinTableHDU, astropy.Table, str
        The table describing the trace_map map. Can be one of: a path to a FITS file
        containing the trace_map infomation, an astropy FITS binary table, or an
        astropy Table object containing the necessary columns

    fits_ext : int
        Default = 1. The FITS extension if a filename for a FITS file is provided
        
    slit_ref : int
        Default = 0. Reference to the slit list index which creates this trace_map. This
        is only used if x is an astropy.Table


    Notes
    -----
    Any table passed to Trace must contain the following columns:
    
    lam : array
        [um] A list of wavelengths which correspond to the coordinate map
        
    x1, [x2,x3], y1 [y2,y3] : array
        [mm] Each an array for the x and y coordinated that trace_map out where
        the spectrum will fall on the detector focal plane. In mm.
        If only x1/y1 are provided, they are considered the central coordinates
        for the slit image. If x2/y2 and x3/y3 are also provided, the slit image is
        projected along the line x1/y1 -> x2/y2 -> x3/y3. This allows effects like
        rotation w.r.t. the trace_map direction and line curvature to be included
    
    weight : array
        [0..1] The weighing of the trace_map. If 100% of the flux is in the trace_map, then
        it is 1. If the flux is split between the trace_maps of two orders, then the relative
        fraction is given here. The array can also contain a wavelength-dependent flux weighting
    
    length : float, optional
        [mm] If x1/y1 and x3/y3 aren't given, this determines the length of the 
        trace on the focal plane.
        
    rotation : float, optional
        [deg] If x1/y1 and x3/y3 aren't given, this determines the rotation of 
        the slit images on the focal plane.
        
        
       
    If a FITS object is passed (filename, astropy.BinTableHDU), the header must contain the 
    following keywords:
    
    CUNIT1, CUNIT2 : [mm] 
        Dimension of spatial coordinates on the focal plane 
        
    CUNIT3 : [um] 
        Dimension of wavelength coordinates in trace_map map
        
    SLIT_REF : int
        Default=0. References the slit the that creates this trace_map
       

    """

    def __init__(self, x, fits_ext=1, slit_ref=0, length=None, rotation=None):
        
        self.params = {"fits_ext" : fits_ext,
                       "slit_ref" : slit_ref,
                       "length"   : length,
                       "rotation" : rotation,
                       "filename" : None}
        
        # holds one/three points which define the slit projection map for a series of wavelengths
        # lam, (x1,y1), (x2,y2), (x3,y3)
        # the units are in millimeters and micron
        
        if isinstance(x, str) and os.path.exists(x):
            self.params["type"] = "FITSFile"
            self.params["filename"] = x
            self.trace_map = Table.read(x, hdu=fits_ext, format='fits')  
            self.slit_ref = self.trace_map.meta["SLIT_REF"]
        
        elif isinstance(x, BinTableHDU):
            self.params["type"] = "BinTableHDU"
            self.trace_map = Table(x.data)
            self.slit_ref = x.header["SLIT_REF"]
            
        elif isinstance(x, Table):
            self.params["type"] = "AstropyTable"
            self.trace_map = x
            self.slit_ref = x.meta["SLIT_REF"] if "SLIT_REF" in x.meta else slit_ref
            
        else: 
            raise ValueError("Need to pass: filename, BinTableHDU, or astropy Table. Not "+type(x))

        # Test if "lam", "x1", "y1", "weight" are in the 
        if not np.all([name in self.trace_map.colnames 
                                    for name in ["lam", "x1", "y1", "weight"]]):
            print(self.trace_map.colnames)
            raise ValueError("The following keywords must be in the Table: lam, x1, y1, weight")


        # ----------------------------------------------------------------------------------
        # Do the funky stuff. Check if x2/y2 and x3/y3 are also there
        if np.all([name in self.trace_map.colnames for name in ["x2", "y2", "x3", "y3"]]):
            
            # For now lets just do a "basic" complex thing where the slit is 
            # only rotated by a constant angle
            rotation_arr = np.arctan2((self.trace_map["y3"]-self.trace_map["y1"]), 
                                      (self.trace_map["x3"]-self.trace_map["x1"]))
            rotation_arr = np.rad2deg(rotation_arr)
            rotation_arr.name = "rotation"
            
            image_length = np.sqrt((self.trace_map["y3"]-self.trace_map["y1"])**2 + \
                                   (self.trace_map["x3"]-self.trace_map["x1"])**2)
            image_length.name = "image_length"
        
            self.trace_map.add_columns((rotation_arr, image_length))
        
        elif np.isscalar(length) and np.isscalar(rotation):
            
            rotation_arr = Column(data=np.array([rotation]*len(self.trace_map)), name="rotation")
            image_length = Column(data=np.array([length]*len(self.trace_map)),   name="image_length")
            
            self.trace_map.add_columns((rotation_arr, image_length))
            
            self.trace_map["x1"].name = "x2"
            self.trace_map["y1"].name = "y2"
            
            dx = 0.5 * length * np.cos(np.deg2rad(rotation))
            dy = 0.5 * length * np.sin(np.deg2rad(rotation))
            
            x1 = self.trace_map["x2"] - dx
            y1 = self.trace_map["y2"] - dy
            x3 = self.trace_map["x2"] + dx
            y3 = self.trace_map["y2"] + dy
            
            x1.name = "x1"
            y1.name = "y1"
            x3.name = "x3"
            y3.name = "y3"
        
            self.trace_map.add_columns((x1, y1, x3, y3, rotation_arr, image_length))
        
        else:
            raise ValueError("Parameters missing. Either [x1/y1, x2/y2, x3/y3] or [x1/y1, length, rotation]")
        
        self.x = self.trace_map["x2"]
        self.y = self.trace_map["y2"]
        self.lam = self.trace_map["lam"]
        
        
    def get_rotation(self, lam):
        """
        Interpolate the rotations array for the required lambda

        Parameters
        ----------
        lam : float
            [um]

        """

        if lam > np.max(self.lam) or lam < np.min(self.lam):
            print(lam, np.min(self.lam), np.max(self.lam))
            raise ValueError("lam outside trace wavelength range")

        rotation = interp2(lam, self.lam, self.trace_map["rotation"])

        return rotation

    def get_weight(self, lam):
        """
        Interpolate the order weights array for the required lambda

        Parameters
        ----------
        lam : float
            [um]

        """

        if lam > np.max(self.lam) or lam < np.min(self.lam):
            print(lam, np.min(self.lam), np.max(self.lam))
            raise ValueError("lam outside trace wavelength range")

        weight = interp2(lam, self.lam, self.trace_map["weight"])

        return weight


    def get_xy(self, lam):
        """
        Interpolate the x,y coordinate arrays for the required line coords

        Parameters
        ----------
        lam : float
            [um]

        Returns
        -------
        x,y : float
            if only one set of coordinates are provided. 

            Or If three coords are available:
        (x,y), (x1,y1), (x3,y3) : float
            centre coords first, then the "left" and "right" coords

        """

        if lam > np.max(self.lam) or lam < np.min(self.lam):
            print(lam, np.min(self.lam), np.max(self.lam))
            raise ValueError("lam outside trace wavelength range")

        x1 = interp2(lam, self.lam, self.trace_map["x1"])
        y1 = interp2(lam, self.lam, self.trace_map["y1"])    

        x2 = interp2(lam, self.lam, self.trace_map["x2"])
        y2 = interp2(lam, self.lam, self.trace_map["y2"])

        x3 = interp2(lam, self.lam, self.trace_map["x3"])
        y3 = interp2(lam, self.lam, self.trace_map["y3"])    

        return [x1, y1], [x2, y2], [x3, y3]


    def is_wavelength_on_chip(self, lam, chip):
        """
        Check if any part of the slit line at a given wavelength inside the given chip area

        Parameters
        ----------
        lam : float
        chip : simcado.detector.Chip

        Returns
        -------
        on_chip : bool
            True if any part of the slit image for this wavelength falls on a chip

        """

        [x1, y1], [x2, y2], [x3, y3] = self.get_xy(lam)

        if np.all(x < chip.x_min for x in [x1,x2,x3]) and \
           np.all(x > chip.x_max for x in [x1,x2,x3]) and \
           np.all(y < chip.y_min for x in [y1,y2,y3]) and \
           np.all(y < chip.y_max for x in [y1,y2,y3]):
            on_chip = False
        else:
            on_chip = True

        return on_chip


    def pixelise_wavelengths(self, pix_size, axis="y", return_xy=False):
        """
        Determine the wavelengths that correspond to a single pixel shift along a given axis

        Parameters
        ----------
        pix_size : float
            [mm] The pixel size on the detector plane.

        axis : str, optional
            Default: "y". The axis along which the discretisation should occur. For vertical traces
            use "y", for horizontal traces give "x".

        return_xy : bool, optinal
            Default is False. If True the x2/y2 coords matching the wavelengths are also returned

        Returns
        -------
        lam_new : array
            [um]

        x_new, y_new : array, optional
            [mm]

        """

        if axis == "y":
            y_new   = np.arange(np.min(self.y), np.max(self.y), pix_size)
            y_new  -= y_new[0] % pix_size    # remove any offset so that y_new is on the pixel centres
            x_new   = interp2(y_new, self.y, self.x)
            lam_new = interp2(y_new, self.y, self.lam)

        if axis == "x":
            x_new   = np.arange(np.min(self.x), np.max(self.x), pix_size)
            x_new  -= x_new[0] % pix_size
            y_new   = interp2(x_new, self.x, self.y)
            lam_new = interp2(x_new, self.x, self.lam)

        if return_xy: 
            return lam_new, x_new, y_new
        else: 
            return lam_new


    def get_wavelength(self, x, y, nearest=True):

        # for a set of focal plane coordinates 

        # are the points inside a trace? point in polygon
        # what is the nearest rotation?
        # where does a line along the rotation direction intersect with the centre line?
        lam = None

        return lam


    def __getitem__(self, key):
        return self.trace_map[key]




class Slit(object):
    """
    Models a Slit
    
    Summary
    -------
    A slit mask  of any dimensions for the first focal plane. The sky will be 
    imaged onto the slit and the slit image is then proejected onto the corresponding
    trace position on the detector focal plane
    
    For a long slit spectrograph only one of these py_objects is needs to be defined.
    The an IFU a slit should be defined for each spaxel
    For a MOS a slit should be definied for each fibre
    
    The is no angle as that is all handles by the spectrograph
    
    Parameters
    ----------
    x, y : float
        [arcsec] The centre of the slit w.r.t the centre of the first focal plane
        
    length, width : float
        [arcsec] All slits are horizontal. Rotation w.r.t to field is handled by 
        the Spectrograph object. 
    
    
    """
    
    
    def __init__(self, x, y, length, width, pix_res=0.004, mask=None, id=None):

        # position in the first focal plane
        # all in [arcsec], with (0,0) being the centre of the plane
        self.x_orig = np.copy(x)
        self.y_orig = np.copy(y)
        
        self.x = x
        self.y = y
        self.length = length   # width
        self.width = width   # height
        #self.angle

        self.pix_res = pix_res
        self.mask = mask if mask.lower() != "none" else None
        self.id = id
    
        # Offsets from atmosphere and anything else that causes a shift
        self.x_offset = 0
        self.y_offset = 0

        # chip object
        self.chip = Chip(x_cen=self.x, y_cen=self.y, 
                         x_len=int(self.length / self.pix_res), 
                         y_len=int(self.width  / self.pix_res),
                         pix_res=self.pix_res, pixsize=0.015,
                         flat_field=self.mask)
        
        
    def apply_offset(self, x_offset=0, y_offset=0, use_orig=True):

        # add an offset and re-make the chip?
        # maybe can just shift the chip coords
        # think of preloading the readout noise

        self.x_offset = x_offset
        self.y_offset = y_offset

        if use_orig:
            self.x = self.x_orig + x_offset
            self.y = self.y_orig + y_offset
        else:
            self.x += x_offset
            self.y += y_offset
        
        ##########################################################
        # I don't like this - remaking the chip is a waste of time
        ##########################################################
        self.chip = Chip(x_cen=self.x, y_cen=self.y, 
                         x_len=int(self.length / self.pix_res), 
                         y_len=int(self.width  / self.pix_res),
                         pix_res=self.pix_res, pixsize=0.015, 
                         flat_field=self.mask)
        
        
#def slits_from_file(filename):


class Spectrograph(OpticalTrain):
    """
    A Spectrograph class, built on top of an OpticalTrain object
    
    Notes
    -----
    The slit knows nothing of the parallactic angle or the field rotation
    This needs to be implemented by rotating the source object 
    offset_tbl holds the offsets from the differential atmospheric dispersion
    
    """
    
    
    def __init__(self, cmds, **kwargs):

        # load normal optical train
        super(Spectrograph, self).__init__(cmds, **kwargs)
        #self.cmds = cmds
        
        # load list of trace_maps
        self.traces = self.get_traces(self.cmds["SPEC_TRACE_LIST"])
        self.n_traces = len(self.traces)
        self.trace_to_slit_ref = [trc.slit_ref for trc in self.traces]

        # load the dimensions of the slits
        self.slits = self.get_slits(self.cmds["SPEC_SLIT_LIST"])
        self.n_slits = len(self.slits)
        self.slit_to_trace_ref = { slit.id : np.where(self.trace_to_slit_ref == slit.id)[0] for slit in self.slits}

        # Get a profile of the offsets induced by the differential atmospheric refraction
        # everything to do with field angle and stuff goes here
        # the slit knows nothing of the angle, but the spectrograph should
        self.offset_tbl = Table()
        self.set_atmo_refraction_offsets(self.cmds["SPEC_REF_WAVELENGTH"], self.cmds["OBS_ZENITH_DIST"])
        
    
        
    def get_slits(self, filename):
        """
        
        """
        
        slit_tbl = ascii.read(filename)
        slit_list = [Slit(row["x_cen"], row["y_cen"], row["length"], row["width"], 
                          row["pix_res"], row["mask"], row["id"]) for row in slit_tbl]
        
        return slit_list


    def get_traces(self, filename):
        """
        import a list of trace_maps from the trace_map_list.tbl file
        deal with the various format that the files may have
        trace_map_list is just a list of Trace py_objects
        
        """
        
        n_traces = fits.getheader(filename, ext=0)["NTRACES"]
        traces   = [Trace(filename, fits_ext=i+1) for i in range(n_traces)]
         
        return traces


    def set_atmo_refraction_offsets(self, lam_ref=None, zenith_dist=None):
        """
        Generates a list of x and y shifts in arcsec due to differential atmospheric diffraction
                
        """

        atmo_offset_lam = np.linspace(self.cmds["SPEC_LAM_MIN"], self.cmds["SPEC_LAM_MAX"], 100)
        if np.isscalar(lam_ref):
            self.cmds["SPEC_REF_WAVELENGTH"] = lam_ref 
        if np.isscalar(zenith_dist) and abs(zenith_dist) < 90: 
            self.cmds["OBS_ZENITH_DIST"] = zenith_dist 
        
        
        effectiveness = self.cmds["INST_ADC_PERFORMANCE"] / 100.

        
        ## get the angle shift for each slice
        lam_shifts = [atmospheric_refraction(lam,
                                                   self.cmds["OBS_ZENITH_DIST"],
                                                   self.cmds["ATMO_TEMPERATURE"],
                                                   self.cmds["ATMO_REL_HUMIDITY"],
                                                   self.cmds["ATMO_PRESSURE"],
                                                   self.cmds["SCOPE_LATITUDE"],
                                                   self.cmds["SCOPE_ALTITUDE"])
                                                   for lam in atmo_offset_lam]
        
        ref_shift = atmospheric_refraction(self.cmds["SPEC_REF_WAVELENGTH"],
                                                 self.cmds["OBS_ZENITH_DIST"],
                                                 self.cmds["ATMO_TEMPERATURE"],
                                                 self.cmds["ATMO_REL_HUMIDITY"],
                                                 self.cmds["ATMO_PRESSURE"],
                                                 self.cmds["SCOPE_LATITUDE"],
                                                 self.cmds["SCOPE_ALTITUDE"])
        
        rel_shift = np.array(lam_shifts) - ref_shift
               
        if np.any(np.abs(rel_shift[-1] - rel_shift[0]) > 10):
            raise ValueError('Diff. Atmo. Refr. too great (>10"), check units')
            
        # Rotate by the parallactic angle
        angle = np.deg2rad(self.cmds["SPEC_SLIT_ANGLE"])
        x = -rel_shift * np.sin(angle) * (1. - effectiveness)
        y = -rel_shift * np.cos(angle) * (1. - effectiveness)

        # return values are in [arcsec]
        self.offset_tbl = Table(data =(atmo_offset_lam, rel_shift, x, y),
                                names=("lam", "abs_shift", "x_offset", "y_offset"))
        
        self.offset_tbl.meta["SPEC_REF_WAVELENGTH"]  = self.cmds["SPEC_REF_WAVELENGTH"]
        self.offset_tbl.meta["SPEC_SLIT_ANGLE"]      = self.cmds["SPEC_SLIT_ANGLE"]
        self.offset_tbl.meta["OBS_ZENITH_DIST"]      = self.cmds["OBS_ZENITH_DIST"]
        self.offset_tbl.meta["INST_ADC_PERFORMANCE"] = effectiveness
        
        
    def get_offsets(self, lam):
        """
        Return the x,y offsets for a given wavelength due to differential atmospheric diffraction
        
        See also
        --------
        set_atmo_refraction_offsets
        
        """
        
        x_offset = np.interp(lam, self.offset_tbl["lam"], self.offset_tbl["x_offset"])
        y_offset = np.interp(lam, self.offset_tbl["lam"], self.offset_tbl["y_offset"])
                             
        return x_offset, y_offset
        




def <Source>.apply_disperser(spec, fpa):

    # 1 find the wavelengths that need to be splattered
    # 2 loop through all wavelengths
        # 3 find the offset of the slit due to atmosphere
        # 4 make a chip that resembles the slit
        # 5 make an image of the slit for the needed wavelength bins
        # 6 find the position on the focal plane where the slit image should be placed
        # 7 find the chips which match this wavelength






def make_zemax_trace_cube(dir_name, return_fits=True, output_name="SPEC_traces.fits"):
    """
    Make a BinTableHDU with the trace maps from ZEMAX files

    """

    pri_hdu = fits.PrimaryHDU()
    pri_hdu.header["NTRACES"] = len(names)
    trace_hdus = [pri_hdu]

    names = [i.split("\\")[-1][:-4] for i in glob.glob(dname+"*.TXT")]

    for name in names:
        fname = name+".TXT"
        order_tbl = read_spec_order(dname+fname)
        order_tbl.remove_columns(("r80_1", "r80_2", "r80_3"))

        weight   = table.Column(0.5*np.ones(len(order_tbl),  dtype=float), name="weight")
        order_tbl.add_column(weight)

        trace_hdu = fits.table_to_hdu(order_tbl)

        trace_hdu.header["DESCRIP"] = name
        trace_hdu.header["CUNIT1"] = "mm"
        trace_hdu.header["CUNIT2"] = "mm"
        trace_hdu.header["CUNIT3"] = "um"
        trace_hdu.header["FILENAME"] = fname
        trace_hdu.header["SLIT_REF"] = 0

        trace_hdus += [trace_hdu]

    trace_hdus = fits.HDUList(trace_hdus)

    if return_fits:
        return trace_hdus
    else:
        trace_hdus.writeto(output_name, clobber=True)


def read_zemax_spec_order(filename):
    """
    Read spectral order definition from a file

    Parameters
    ----------
    filename : str


    Returns
    -------
    An `astropy.Table` with columns
        - lam : wavelength
        - x_1, x_2, x_3 : Column numbers of left edge, centre and right edge
               of lines of constant wavelength
        - y_1, y_2, y_3 : Row numbers of left edge, centre and right edge
               of lines of constant wavelength
        - r80_1, r80_2, r80_3 : radii of 80% encircled energy

    Notes
    -----
    The orders file has the following format:
    - two lines of information
    - for each wavelength six lines:
      - a line with "index= ... wavelength= ".
      - a flag (0 is good)
      - a line each for the left edge, centre and right edge of the 2D
        trace_map at the given wavelength. Format is "X1= ... Y1= ... r(EE80)= ..."
      - a separator line

    """

    # Read file
    with open(filename) as fp1:
        lines = fp1.readlines()

    nlines = len(lines)

    # Initialize lists
    lam = []
    x_1 = []
    x_2 = []
    x_3 = []
    y_1 = []
    y_2 = []
    y_3 = []
    r80_1 = []
    r80_2 = []
    r80_3 = []

    # Extract the order number
    lhs = lines[1].split(',')[0]
    order = int(float(lhs.split()[2]))

    # Drop the first two lines
    iline = 2

    # Loop over group of six lines at a time
    while iline < nlines - 1:
        # Use only "good" lines
        flag = float(lines[iline + 1])
        if flag == 0:
            # shortcuts for the next four lines
            lamline = lines[iline]
            x1line = lines[iline + 2]
            x2line = lines[iline + 3]
            x3line = lines[iline + 4]

            # Extract information
            lam.append(float(lamline.split()[3]))
            x_1.append(float(x1line.split()[1]))
            y_1.append(float(x1line.split()[3]))
            r80_1.append(float(x1line.split()[5]))
            x_2.append(float(x2line.split()[1]))
            y_2.append(float(x2line.split()[3]))
            r80_2.append(float(x2line.split()[5]))
            x_3.append(float(x3line.split()[1]))
            y_3.append(float(x3line.split()[3]))
            r80_3.append(float(x3line.split()[5]))
        iline += 6

    # return as numpy arrays
    from astropy.table import Table
    order_table = Table([lam, x_1, x_2, x_3, y_1, y_2, y_3, r80_1, r80_2, r80_3],
                        names=['lam', 'x1', 'x2', 'x3', 'y1', 'y2', 'y3',
                               'r80_1', 'r80_2', 'r80_3'],
                        meta={'Order': order})

    return order_table


def interp2(x, xp, fp):
    """
    The same as numpy.interp but accepts decreasing xp values

    See also
    --------
    numpy.interp

    """

    if np.all(xp[1:] - xp[:-1] > 0):
        f = np.interp(x, xp, fp)
    elif np.all(xp[1:] - xp[:-1] < 0):
        f = np.interp(x, xp[::-1], fp[::-1])
    else:
        raise ValueError("xp does not continuously increase or decrease")

    return f

    
    
def overlay_image(small_im, big_im, coords, sub_pixel=False):
    """
    Overlay small_im on top of big_im at the position specified by coords

    ``small_im`` will be centred at ``coords``

    Adapted from:
    ``https://stackoverflow.com/questions/14063070/overlay-a-smaller-image-on-a-larger-image-python-opencv``
    
    """

    # TODO - Add in a catch for sub-pixel shifts
    if sub_pixel:
        raise NotImplementedError
    
    x, y = np.array(coords, dtype=int) - np.array(small_im.shape) // 2
    
    # Image ranges
    x1, x2 = max(0, x), min(big_im.shape[0], x + small_im.shape[0])
    y1, y2 = max(0, y), min(big_im.shape[1], y + small_im.shape[1])
    
    # Overlay ranges
    x1o, x2o = max(0, -x), min(small_im.shape[0], big_im.shape[0] - x)
    y1o, y2o = max(0, -y), min(small_im.shape[1], big_im.shape[1] - y)
    
    # Exit if nothing to do
    if y1 >= y2 or x1 >= x2 or y1o >= y2o or x1o >= x2o:
        return big_im

    big_im[x1:x2, y1:y2] += small_im[x1o:x2o, y1o:y2o]   
    
    return big_im   