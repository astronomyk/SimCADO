###############################################################################
# UserCommands
#
# DESCRIPTION
#
# UserCommands is essentially a dictionary that holds all the variables that
# the user may wish to change, plus some methods to get certain subclasses
# out, i.e. anything that starts with "ATMO"
#
#
# Classes:
#   
#
# Methods:
#   
#
#



#=== All the current variables for playing with ===

##### DIRECTORY STUCTURE #####
HOME_DIR               Code
TRANS_CURVE_DIR        ../Transmission_curves/
EMISS_CURVE_DIR        ../Emission_curves/
PSF_FILE_DIR           ../PSFs/
INPUT_CUBE_DIR         ../Input_data


##### ATMOSPHERE #####
ATMO_USE_SKYCALC_TC     yes         # use a skycalc FITS table, or provide separate tables for tc/ec
ATMO_USE_SKYCALC_EC     yes         # use a skycalc FITS table, or provide separate tables for tc/ec
ATMO_SKYCALC_FILE      ../Emission_curves/skytable.fits
ATMO_TC                none        # transmission curve file path
ATMO_EC                none        # emission curve file path
ATMO_BG_PHOTONS        10        # [ph/s] for the bandpass in use

ATMO_TEMPERATURE       0           # deg Celcius
ATMO_PRESSURE          750         # millibar
ATMO_REL_HUMIDITY      60          # %
ATMO_ZENITH_DIST       60          # deg from zenith



##### TELESCOPE ######
SCOPE_ALTITUDE         3060        # meters above sea level
SCOPE_LATITUDE         -24.589167  # decimal degrees
SCOPE_LONGITUDE        -70.192222  # decimal degrees

SCOPE_NUM_MIRRORS      6           # number of reflecting surfaces
SCOPE_M1_DIAMETER_OUT  39.3        # meters
SCOPE_M1_DIAMETER_IN   4.0         # meters
SCOPE_TC_MIRROR        ../Transmission_curves/TC_mirror_mgf2agal.dat
SCOPE_EC_MIRROR        none

SCOPE_USE_PSF_FILE     no          # import a PSF from a file or generate PSF internally
SCOPE_AO_EFFECTIVENESS 100         # percentage of seeing PSF cancelled out AO - 100% = diff limited, 0% = 0.6" seeing
SCOPE_JITTER_FWHM      5           # assuming a gaussian for general jitter, FWHM in milli-arcsec


##### INSTRUMENT #####
INST_TEMPERATURE       -190        # deg Celsius - inside temp if instrument 
INST_NUM_MIRRORS       8           # number of reflecting surfaces
INST_TC_MIRROR         ../Transmission_curves/TC_mirror_mgf2agal.dat
INST_TC_ENTR_WINDOW    none
INST_TC_DICHROIC       none
INST_TC_FILTER         ../Transmission_curves/TC_fitler_Ks.dat
INST_TC_ADC            none
INST_ADC_EFFICIENCY    0          # %
INST_ADC_NO_SURFACES   4



##### DETECTOR #####
FPA_DARK_MEDIAN        0.01       # e-/s/px
FPA_DARK_STDEV         0.01       # e-/s/px
FPA_READOUT_MEDIAN     5          # e-/px
FPA_READOUT_STDEV      1          # e-/px

FPA_QE                 ../Transmission_curves/TC_detector_H4RG.dat
FPA_DISTORTION_MAP     none
FPA_LINEARITY_CURVE    none
FPA_GAIN               1          # e- to ADU conversion
FPA_WELL_DEPTH         100000     # number of photons collectable before pixel is full


##### OBSERVATIONS #####
OBS_ALT                none
OBS_AZ                 none
OBS_SEEING             none
OBS_EXPTIME            600         # [sec] simulated exposure time

OBS_USE_FILTER_LAM     yes        # yes/no to basing the wavelength range off the filter non-zero range - if no, specify LAM_MIN, LAM_MAX
OBS_LAM_MIN            1.9        # [um] lower wavelength range of observation
OBS_LAM_MAX            2.41       # [um] upper wavelength range of observation
OBS_LAM_BIN_WIDTH      0.1        # [um]    
OBS_PIXEL_SCALE        4          # [mas]
OBS_PIXEL_THRESHOLD    0        # photons per pixel summed over the wavelength range

OBS_INPUT_FILENAME     /media/sf_Share_VW/Data_for_SimCADO/sim_input/cluster10000_1.0_Myr_800.0_kpc.fits
OBS_FITS_EXT		   0          # the extension number where the useful data cube is

OBS_OUT_SIGNAL         yes
OBS_OUT_NOISE          yes
OBS_OUT_READOUT        yes
OBS_OUT_REDUCED        yes