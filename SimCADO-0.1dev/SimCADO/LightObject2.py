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
from astropy.convolution import convolve, convolve_fft
import numpy as np
import scipy.ndimage.interpolation as spi

import os
import warnings

class LightObject(object):

    def __init__(self, **kwargs):
    """
    Keywords:
    - lam: [µm] an array of size L with wavelengths for the spectra
    - x, y: [pix] arrays of size N holding the pixel coordinate information for
            N sources
    - spectra: [photons] a 2D array of size (S,L) with a spectrum for S unique 
                sources
    - spec_ref: [int] an array holding N references joining each of the N source
                pixel coordinates to one of the S unique spectra
    - weights: [float] an array of size N with weights for each source
    
    
    """
    
    
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
        
        self.lam      = np.asarray(kwargs["lam"])
        self.spectra  = np.asarray(kwargs["spec_list"])
        self.x        = np.asarray(kwargs["x"])
        self.y        = np.asarray(kwargs["y"])
        self.spec_ref = np.asarray(kwargs["spec_ref"])
        self.weights  = np.asarray(kwargs["pixel_weights"])
               
        # add a second dimension to self.spectra so that all the 2D calls work
        if len(self.spectra.shape) == 1:
            self.spectra.shape = np.asarray([self.spectra.shape]*2)
        
        self.array = np.zeros((self.params("NAXIS1"),
                               self.params("NAXIS2")), dtype=np.float32)
        
    def __repr__(self):
        return self.info['description']
    
    def __array__(self, x):
        return self.array
    
    def poissonify(self):
        """ Add a realisation of the poisson process to the photon signal """
        self.array = np.random.poisson(self.array)
    
    
    def apply_psf_cube(self, psf_cube, lam_bin_edges, sub_pixel=False, 
                       export_slices=False):
        """
        For all PSFs in a PSFCube, generate an array containing the number of 
        photons expected in the wavelength range of that PSF, and where they
        will land on the FOV. Then apply the PSF to this "ideal" FOV.
        
        ??? Issue
        zoom_res for get_slice_photons is set at 10. Should this be variable?
        
        Keywords:
        - psf_cube
        - lam_bin_edges
        
        Optional keywords:
        - sub_pixel: [False, True] simulate sources that aren't directly in the
                     centre of a pixel
        - export slices: [False, True] export all the slices to a fits file
        """

        if len(lam_bin_edges) != len(psf_cube) + 1:
            warnings.warn("Number of PSFs does not fit to lam_bin_edges")
    
        x_int, y_int = self.x.astype(int), self.y.astype(int)
        slice_photons = np.zeros((len(self.lam))
        slice_array = np.zeros((self.params("NAXIS1"),
                                       self.params("NAXIS2")), dtype=np.float32)
        
        for i in range(len(psf_cube)):
            
            psf = psf_cube[i]
            lam_min, lam_max = lam_bin_edges[i], lam_bin_edges[i+1]
            slice_photons = self.get_slice_photons(lam_min, lam_max, zoom_res=10)
            
            # if sub pixel accuracy is needed, be prepared to wait. For this we
            # need to go through every source spectra in turn, shift the psf by
            # the decimal amount given by pos - int(pos), then place the a 
            # certain slice of the psf on the output array. 
            if sub_pixel:
                dx, dy = self.x - x_int, self.y - y_int
                ax, ay = np.array(slice_array.shape)//2           
                bx, by = np.array(psf.array.shape)//2
                
                for p in range(len(slice_photons)):
                    psf_tmp = spi.shift(psf, (dx[p],dy[p]), order=1)
                    x_pint, y_pint = x_int[p], y_int[p]
                    
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

                    slice_array[ax0:ax1, ay0:ay1] = psf[bx0:bx1, by0:by1] \
                                            * slice_photons[p] * self.weights[p]
                    
            else:
                # If astrometric precision is not that important and everything 
                # has been oversampled, use this section.
                
                slice_array[x_int, y_int] = slice_photons * self.weights
                slice_array = convolve_fft(self.slice_array, psf.array)
            
            self.array += slice_array
            
            if export_slices:
                filename = "../data/tmp.fits"
                if not os.path.exists(filename):
                    hdu = fits.PrimaryHDU(slice_array)
                    hdu.writeto(filename)
                else:
                    hdu = fits.open(filename, mode="append")
                    hdu.append(fits.ImageHDU(slice_array))
                    hdu.update_extend()
                    hdu.flush()
                os.rename("../tmp_data/tmp.fits", "../tmp_data/slices.fits")
    
    
    def apply_plane_effect(self, plane):
    
    
    def apply_spectral_curve(self, spec_curve):
        
        
    
    
    
    
    
    
    def get_slice_photons(self, lam_min, lam_max, zoom_res = 10):
        """
        Caluclate how many photons for each source exist in the wavelength bin
        defined by lam_min and lam_max.
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
        
        #spec_zoom = np.asarray([np.interp(lam_zoom, lam[i0:i1], spec[i0:i1]) for spec in spectra])
        for i in range(len(spectra)): 
            spec_zoom[i,:] = np.interp(lam_zoom, lam[i0:i1], spectra[i,i0:i1])

        slice_photons = np.trapz(spec_zoom, lam_zoom, axis=1)  
        return slice_photons

        