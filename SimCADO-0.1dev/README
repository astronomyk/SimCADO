###############################################################################
#                                 SimCADO                                     #
#                 Kieran Leschinski and Oliver Czoske                         #
#                        Last updates 29.01.2016                              #
###############################################################################

###############################################################################
WHAT IS SIMCADO?

SimCADO is the instrument data simulator for MICADO. It aims to simulate the 
complete optical train from the source of the photons - i.e. the data cube 
provided by the user - through the atmosphere, E-ELT, MICADO, detector array 
and into the 1s and 0s of a readout fits file. 


###############################################################################
PREREQUISITES FOR SIMCADO
-Python 2.7 or >3.4
- Numpy
- Scipy
- Matplotlib
- Astropy


###############################################################################
HOW DO I USE SIMCADO v0.2?

Before getting too excited, please read the section: WHAT CAN SIMCADO v0.2 NOT 
DO? This will save everyone time.

That said, SimCADO is easy to use. It is currently a simple Python2.7 script 
that accepts as a command line argument a single user written configuration 
file. To use SimCADO v0.2, enter the /SimCADO/Code/ folder and type:

$ python  simcado.py  my_config_file.txt

A word of warning. All the default parameters are in the file simcado.default. 
This file should be used as a reference for the user's own config file only.
NEVER change simcado.default, unless you know what you're doing. SImCADO reads 
in all the other configuration data from this file for every run before looking 
at the user's config file. Anything in the user's config file takes precedence 
and the default values are overwritten on the fly.

All the default instrument data is bundled with SIMCADO. If you would like to 
use your own, feel free.


###############################################################################
WHAT HEADER KEYWORDS DOES MY DATACUBE NEED?

SimCADO will try to read any FITS file and take what information it can find 
from the header. However FITS files are like snowflakes - no two are ever the 
same. Please make sure that your data cube has the following header keywords so 
that SimCADO can deliver real results, rather than just a meaningless pile of 
1s and 0s.

NEEDED KEYWORDS    DESCRIPTION
------------------------------------
NAXIS	           3 - needs to be a cube at the moment
NAXISn             < 4096
CRVALn
CRPIXn
CDELTn

HELPFUL KEYWORDS   DESCRIPTION
------------------------------------
EXPTIME            [seconds] exp. time of cube, default 1.s
BUNIT              if no units are given, defaults to [ph/s]
ADU                see GAIN
GAIN               [e-/ADU] the default is 1.
PIX_RES            [mas] pixel resolution, default is 4.mas

All other keywords in the header are ignored at the moment


###############################################################################
GETTING STARTED WITH SIMCADO v0.2

In this demonstration version I have include a file with the most frequently 
used parameters - my_config_file.txt 
Feel free to edit this file, copy it, delete keywords that you don't use or 
copy in parameters from simcado.default.

Below I have a short description of the most used keywords

FREQUENTLY USED KEYWORDS
-------------------------------------------------------------------------------
If yes, use the photon counts from the skycalc model for the sky background 
emission. If no, assume ATMO_BG_PHOTONS as the number of sky background photons 
per second coming through the filter bandpass

ATMO_USE_SKYCALC_EC    yes       # [yes/no]
ATMO_BG_PHOTONS        50        # [ph/s]

Because we don't yet have PSFs from the SCAO and MAORY teams, we have to make 
do with an approximation of the telescope PSF and the effectiveness of the AO.
100% = diffraction limited PSF, 0% = 0.6" seeing PSF. Furthermore, any extra 
gaussian PSF contribution can be added in with SCOPE_JITTER_FWHM

SCOPE_AO_EFFECTIVENESS 100       # [0...100] percentage
SCOPE_JITTER_FWHM      1         # [mas]

Inside MICADO we can specify the filter to use. Look in the folder to see which 
filters are available, or use your own. We also don't have a proper model of 
the ADC, so here we use a basic efficiency percentage. 

INST_TC_FILTER         ../Transmission_curves/TC_filter_J.dat
INST_ADC_EFFICIENCY    80         # [0...100] percentage

Basic detector noise characteristics, such as readout noise, the gain of the 
detector and the well depth of the pixels.

FPA_READOUT_MEDIAN     5          # [e-/px]
FPA_GAIN               1          # e- to ADU conversion
FPA_WELL_DEPTH         100000     # [1...1E99]

OBS_EXPTIME is time SimCADO should "observe" the source cube. The 
OBS_INPUT_FILENAME points to the 3D FITS cube that will be "observed". 
OBS_FITS_EXT is the extension number for the useful part of the cube. (There's 
no point to reading in a FITS table in extension 0!)

OBS_EXPTIME            300        # [sec]
OBS_INPUT_FILENAME     ../Input_data/My_Amazing_Cube.fits
OBS_FITS_EXT		   0          # [0...n]

As the input cubes should be noiseless, yet SimCADO will accept real data (i.e. 
KMOS or MUSE cubes) the background noise should not be taken as "signal". 
Therefore we set OBS_PIXEL_THRESHOLD to just above the noise level in the image 
so that only the "real" source photons are used as input for SimCADO. If set to 
0, all the original image noise will added to the detector array output image.

OBS_PIXEL_THRESHOLD    100        # [0...n]

OBS_LAM_BIN_WIDTH is the spectral resolution for the simulation. This should be 
worse (by a factor of >2) than the spectral resolution of the input cube. Long 
story. In essence using 0.1um for JHK if the ADC is at 100% is perfectly fine. 
If the ADC is at 0%, spectral resolution in J-band should be ~0.025um to avoid 
star splotches.

OBS_LAM_BIN_WIDTH      0.1        # [um]    

There are 4 possible output files. Each takes the file name
of the input cube plus a suffix before the .fits extension:
SIGNAL -  is an image of the pure source cube signal arriving at the detector 
          (*.signal.fits)
NOISE -   is the noise generated by the detector and the sky background 
          emission (*.noise.fits)
READOUT - is the combination of SIGNAL and NOISE as seen by the detector array 
		  (*.fpa.fits)
REDUCED - is READOUT minus the median background flux (*.red.fits)
		 
OBS_OUT_SIGNAL         no         # [yes/no]
OBS_OUT_NOISE          no         # [yes/no]
OBS_OUT_READOUT        yes        # [yes/no]
OBS_OUT_REDUCED        no         # [yes/no]

###############################################################################
WHAT CAN SIMCADO v0.2 DO?

Big note of warning - SimCADO v0.2 is not yet finished. When the full version 
is available, you will know. We will be making a lot of noise about it. In the 
meantime we offer SimCADO v0.2 as a demonstration of what we are currently 
working towards. 

That said, SimCADO v0.2 is able to simulate realistic looking detector array 
FITS images for various configurations of the optical train. Features include:
- IzJHK filters + user defined filters
- basic ADC 
- basic AO
- detector array, pure signal, and pure noise output images
- variable instrument mirror configuration
- variable detector noise, background noise 
...

SimCADO currently needs a 3D data cube to generate detector readout images. 
Functionality for working with 2D images will be implemented very soon.

###############################################################################
WHAT CAN SIMCADO v0.2 NOT DO?

NOTE:
Currently ONLY imaging mode is offered. Once the Imaging mode is bug free, we 
will begin implementing the others.

NOTE: 
SimCADO v0.2 only works on Python 2.7. We apologise that we have not yet made 
the leap into Python 3. This will come with the SimCADO v1.0 release.

As SimCADO is still in its pre-alpha stage, there are several things it can't 
do yet. All of these points are on our to-do list, so please be patient. 
- read 2D images
- use SCAO or MCAO PSF files
- distort the image (sky rotation, instrument distortion, )
- model the true detector characteristics
- use oversampled data


###############################################################################
EXPLANATION OF ALL THE KEYWORDS

Coming soon to a laptop near you!