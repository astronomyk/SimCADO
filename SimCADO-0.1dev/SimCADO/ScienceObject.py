###############################################################################
# ScienceObject
#
# DESCRIPTION
# The ScienceObject should appear in the form of a FITS cube which stays on
# disk and is never fully read into RAM. A 4x oversampled 4k x 4k image with
# using float32 in python takes up 2GB of RAM. Therefore the simplest way to 
# avoid memory problems is to never read in the full data cube to begin with.
# We can get around this problem by reading in the spectral channels one by one.
# The ScienceObject must therefore only contain a pointer to a copy of the 
# relevant layers of the FITS file, which can be manipulated and collapsed.
#
# ScienceObject will need to contain spatial and wavelength information:
# - pix_res
# - pix_width
# - pix_height
# - lam_res
# - lam_bin_centers
# - lam_bin_edges
#
# ScienceObject will be in photons/s/m2/pixel_fov/wavelength_slice. Therefore
# ScienceObject will also need to know the telescope collecting area and 
# integration time to convert the cube into photons/spaxel. The following info
# should be pulled from the UserCommands dictionary:
# - area
# - exptime
#
#
#
#
#
#
#
#
# spectral layers 
#
#
#
#
#
#
# Classes:
#  LightCube(object)
#  
#
# Methods:
#   
#
#
#
#
#



from astropy.io import fits

## These classes and functions are exported to the package
__all__ = ["LightCube"]

class LightCube(object):
    """Cube holding the light passing through the optical system"""

    def __init__(self, fname, extnum=0):
        self.info = dict([])
        self.info['Filename'] = os.path.basename(fname)
        self.info['Path'] = os.path.dirname(fname)
        
        fp = fits.open(fname)
        self.header = fp[extnum].header
        fp.close()
        

    def __repr__(self):
        string = "LightCube\nFile name: %s" % self.info['Filename']
        return string
    
    def apply_atmosphere(self, atmosphere_model):
        self.info['atmosphere'] = 'applied'

    def apply_psf(self, psf_cube):
        self.info['psf'] = 'applied'

    def apply_throughput(self, throughput):
        self.info['throughput'] = 'applied'

    def collapse(self):
        self.info['collapse'] = 'applied'

    def apply_geometry(self):
        self.info['geometry'] = 'applied'

    def make_detector_image(self, detector):
        self.info['detector_image'] = 'applied'

    def add_background(self, bg_level):
        self.info['background'] = 'applied'

        
