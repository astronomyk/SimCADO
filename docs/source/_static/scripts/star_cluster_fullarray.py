#!/usr/bin/env python

"""
We create a star cluster that mostly fills all 9 detectors. We use a 
MAORY-like J-band PSF that is uniform over the whole filed. We set the exposure time so
the brightest sources have a peak of about 10000-20000 counts. 

The PSF used here is a MAORY-like PSF with a Strehl of 7% (PSF_MCAO_J_Strehl7.fits). 
It must be located where SimCADO can find it
"""

import numpy as np
from astropy.io import fits
import simcado
print(simcado.__data_dir__)


#  Generating the star cluster
#  simcado.source.cluster generates a star cluster with age=0. As such it produces a few very 
# bright stars that will peak at high photon counts. To have more visible sources I added few of 
# these clusters thus generating a stellar cluster with a mass of 1.5e5 Msum

src = simcado.source.cluster(mass=1E4, distance=2e6,  half_light_radius=200)
for i in range(15):
    src = src + simcado.source.cluster(mass=1E4, distance=3e6,  half_light_radius=250)

t_exp = 10  #s
sim_hdu = simcado.run(src, OBS_EXPTIME=t_exp, detector_layout="full",
                      INST_FILTER_TC="TC_filter_J.dat", 
                      SCOPE_PSF_FILE="PSF_MCAO_J_Strehl7.fits",
                      FPA_LINEARITY_CURVE=None)


for i in range(9):
    print("Maximum value of chip", i+1, ":",  np.max(sim_hdu[i+1].data))

# Saving the file
filename = "star_cluster.fits"
print("Saving the results to", filename)
sim_hdu.writeto(filename, overwrite=True)
