# 1. FVPSF should return the PSF for a position in the FOV and a given lambda
# 2. should throw errors when:
#   - file doesn't exist
#   - file doesn't have a image or table in the 1st extension
# 3. should have attributes:
#   - lam_bin_centers : pulled from the header of the hduNs
#   - layer_map : an ImageHDU with a map of where each layer is valid
#   - layer_table : The BinTableHDU if available
#   - _make_layer_map : makes a layer_map from a BinTableHDU is needed
#   - mask(wave, pos) : returns an ImageHDU with WCS with a mask of the valid
#                       region for a given wavelength and position in the FOV
#   - shape : (N_EXT, N_LAYER, NAXIS2, NAXIS1)
#   - nearest(wave, pos=None, hdu=False) : should return the array for the
#                                          given position and wavelength
#   - defaults : a dictionary with {"wave": , pos: (,)} so that .array can be
#                used for backwards compatibility
#   - set_defaults(wave, pos)
#   - array : returns an array for the defaults values of wave and pos
#   - psf : returns self, as array returns a layer based on defaults

