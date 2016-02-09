###############################################################################
# PlaneEffect
#
# DESCRIPTION
# The PlaneEffect object is used to simulate effects that occur on all spectral
# layers equally, for example, sky rotation, telescope jitter, distortion, etc.
# To do this in the most general way possible, a PlaneEffect object contains 
# 3 planes representing the deviation in position from an ideal optical train.
# The values in each of the 3 planes represent the distance the pixel should 
# move in the x and y directions, and a weighting value.
#
# Several subclasses generate the various effects that occur. For example: 
# - Rotation
# - Distortion
# - Translation
# - FlatField
#
#
#
# Classes:
#  PlaneEffect
#
# Subclasses:
#  Rotation(PlaneEffect)
#  Distortion(PlaneEffect)
#  Translation(PlaneEffect)
#  FlatField(PlaneEffect)
#
# Methods:
#
#

class PlaneEffect(object):
    """
    
    """

    def __init__(self, **kwargs):
    
        self.info = dict([])
        self.info['created'] = 'yes'
        self.info['description'] = "Point spread function (single layer)"
    
        self.x_map = np.zeros((2,2))
        self.y_map = np.zeros((2,2))
        
    

class 
