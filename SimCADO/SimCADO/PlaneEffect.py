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
