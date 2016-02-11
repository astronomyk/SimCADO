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
import scipy.ndimage.interpolation as spi


class LightObject(object):

    def __init__(self, **kwargs):

    # - lam_res
    # - lam_bin_centers
    # - lam_bin_edges
        
        self.params = { "pix_res"   :0.004,
                        "NAXIS1"    :4096,
                        "NAXIS2"    :4096,
                      }
        self.params.update(kwargs)

        self.info = dict([])
        self.info['created'] = 'yes'
        self.info['description'] = "List of spectra and their positions"
        
        self.lam       = np.asarray(kwargs["lam"])
        self.spectra   = np.asarray(kwargs["spec_list"]))
        if len(self.spectra.shape) == 1:
            self.spectra.shape = np.asarray([self.spectra.shape]*2)
        
        self.positions = np.asarray((kwargs["x"], kwargs["y"], 
                                kwargs["spec_ref"], kwargs["pixel_weights"]))
        
        self.working_slice = np.zeros((self.params("NAXIS1"),
                                       self.params("NAXIS2")), dtype=np.float32)
        self.array = np.zeros((self.params("NAXIS1"),
                               self.params("NAXIS2")), dtype=np.float32)
        
    def __repr__(self):
        return self.info['description']
    
    def __array__(self, x):
        return self.array
    
    def get_slice_photons(self, lam_min, lam_max, zoom_res = 10):
        """
        
        """
        
        # Check if the slice limits are within the spectrum wavelength range
        if lam_min > self.lam[-1] or lam_max < self.lam[0]:
            print((lam_min, lam_max), (self.lam[0], self.lam[-1]))
            raise ValueError("lam_min or lam_max outside wavelength range of spectra")

        # find the closest indices i0, i1 which match the limits 
        x0, x1 = np.abs(self.lam - lam_min), np.abs(self.lam - lam_max)
        i0, i1 = np.where(x0 == np.min(x0))[0][0], np.where(x1 == np.min(x1))[0][0]
        if self.lam[i0] > lam_min and i0 > 0: 
            i0 -= 1 
        if self.lam[i1] < lam_max and i1 < len(self.lam): 
            i1 += 1 

        zoom_factor = zoom_res * (i1 - i0)
        lam_zoom  = np.linspace(lam_min, lam_max, zoom_factor)
        spec_zoom = np.zeros((spectra.shape[0], zoom_factor))
        for i in range(len(spectra)): 
            spec_zoom[i,:] = np.interp(lam_zoom, lam[i0:i1], spectra[i,i0:i1])

        photons = np.trapz(spec_zoom, lam_zoom, axis=1)  
        return photons
    
    def apply_spectral_curve(self, spec_curve):
     
        
    def convolve(self, psf):
  
    