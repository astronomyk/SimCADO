###############################################################################
# SimulationRun
#
# Deliverables
#
#
# Need 2 different UserCommands dictionaries as input
#  - Observation parameters (Alt, Az, Exptime, etc)
#  - Observatory parameters (area, instrument configuration, psfs, etc)
#
#
#
# Classes:
#
#
# Methods:
#
#
# 1. Read in UserCommands
#    Create UserCommands objects for: optical train and observing run
#    The observing run will be used for science, atmo and mirror photons
#
# 2. Create OpticalTrain, 
#    which in turn creates all the individual objects
#
# 3. Create 
