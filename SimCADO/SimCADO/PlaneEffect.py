###############################################################################
# PlaneEffect
#
# DESCRIPTION
# The PlaneEffect functions are used to simulate effects that occur on all spectral
# layers equally, for example, sky rotation, telescope jitter, distortion, etc.
# To do this in the most general way possible, a PlaneEffect function contains 
# 3 planes representing the deviation in position from an ideal optical train.
# The values in each of the 3 planes represent the distance the pixel should 
# move in the x and y directions, and a weighting value.
#
# Several functions generate the various effects that occur. For example: 
# - Rotation
# - Distortion
# - Translation
# - FlatField
# Some PlaneEffects only need to act on the positions of the incoming photons,  
# e.g. ADC, while others are applicable to the whole array, e.g. Distortion, 
# Flat field. As each PlaneEffect
# - CoordEffect
# - ArrayEffect
# 
# As each PlaneEffect 
#
# Classes:
#  CoordEffect
#  ArrayEffect
#
# Subclasses:
#  Rotation(ArrayEffect)
#  Distortion(ArrayEffect)
#  FlatField(ArrayEffect)
#  Translation(CoordEffect)

#
# Methods:
#
#
#
from copy import deepcopy

import numpy as np
import scipy.ndimage as spi


class CoordEffect(object):
    """
    """

    def __init__(self, dx, dy):
        self.dx = dx
        self.dy = dy

    def apply(self, x, y)
        return x + self.dx, y + self.dy
        
        
        
class ADC_Effect(CoordEffect):
    
    def __init__(lam, angle = 0, **kwargs):
    
        shift = atmospheric_refraction(lam, **kwargs):
        dx = shift * np.cos(angle / 57.29578)
        dy = shift * np.cos(angle / 57.29578)

        super(ADC_Effect, self).__init__(dx, dy)

class ArrayEffect(object):
    
    def __init__(self, x, y, weight):
        self.x = np.zeros()
    
    
    def apply(self, array)
    
    


        
        
def rotate_array(x, y, angle, center):
    pass

    

    return x+dx, y+dy

def shift_array(arr, dx, dy):
    pass

def distort_coords(x, y, x_arr, y_arr)
    pass

def distort_array(arr, x_arr, y_arr, weight_arr)
    pass


        
    def rotate(x, y, angle, center):
        pass
    
    def shift(dx, dy):
        pass

















