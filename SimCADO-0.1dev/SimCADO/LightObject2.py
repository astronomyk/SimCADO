###############################################################################
# LightObject
#
# DESCRIPTION
# The LightObject is essentially a list of spectra and a list of positions. The
# list of positions contains a reference to the relevant spectra. The advantage
# here is that if there are repeated spectra in a data cube, we can reduce the 
# amount of calculations needed. Furthermore, if the input is originally a list
# of stars, etc where the position of a star is not always and integer multiple 
# of the plate scale, we can keep the information until the PSFs are needed.
#
# The LightObject contains two arrays:
#  - PositionArray: 
#  - SpectrumArray



# Flow of events
# - Generate the lists of spectra and positions
# - Apply the transmission curves [SpectrumArray]
# - shrink the 1D spectra to the resolution of the PSFCube layers [SpectrumArray]
# - apply any 2D PlaneEffects [PositionArray]
# for i in len(slices)
#   - generate a working slice [PositionArray, SpectrumArray, WorkingSlice]
#   - Apply the PSF for the appropriate wavelength [WorkingSlice]
#   - Apply any wavelength dependent PlaneEffects [WorkingSlice]
#   - apply Poisson noise to the photons in the slice [WorkingSlice]
#   - add the WorkingSlice to the FPA [WorkingSlice, FPArray]



from astropy.io import fits
import numpy as np


class LightObject(object):

    def __init__(self, **kwargs):

        self.spectra   = SpectrumArray(kwargs["lam"], kwargs["spec_list"])    
        self.positions = PositionArray(kwargs["x"], kwargs["y"], 
                                        kwargs["spec_ref"], kwargs["weight"])

        

        
    def generate_slice(self, lam_min, lam_max)
    
        return slice
    
    
class SpectrumArray(object):

    def __init__(self, lam, spec_list):
        self.lam = lam
        self.array = np.asarray(spec_list)
        
    
    def shrink_spectra(self, lam_bin_edges, action="sum"):

        return low_res_spectra
    
    
class PositionArray(object):

    def __init__(self, x, y, spec_ref, weight):

    
    

class WorkingSlice(object):

    def __init__(self, spec_array, pos_array, lam_bin_edges, pix_res)
    
    
    
    
    
    
    
    
    

    
    