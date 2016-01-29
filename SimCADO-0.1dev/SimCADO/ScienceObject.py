###############################################################################
# ScienceObject
#
# DESCRIPTION
#
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

        
