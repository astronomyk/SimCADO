#!/usr/bin/env python

"""
I simulate 20 stars between H=20-24 magnitudes (random) with positions defined by a normal distribution with FHWM=1". The central position of the distribution changes and it is moved on the Y-axis by 7, 14 and 21 arcsec. So, we can see the PSF change when going off axis. Single exposure in the H-band with Texp=

This experiment needs AnisoCADO PSF. They can be downloaded from https://simcado.readthedocs.io/en/latest/user_docs/9_PSFs.html  and should be made available to SimCADO by e.g. placing them in the SIM_DATA_DIR
"""

import numpy as np
import simcado
print(simcado.__data_dir__)


# Creating the sources
magnitudes = np.random.uniform(low=20, high=24, size=20) # Vega
x_pos = np.random.normal(scale=1., size=20) 
y_pos = np.random.normal(scale=1., size=20)  


point_src1 = simcado.source.stars(mags=magnitudes, x=x_pos, y=y_pos, 
                                     filter_name='TC_filter_H.dat', spectype="A0V")

point_src2 = simcado.source.stars(mags=magnitudes, x=x_pos, y=y_pos+7, 
                                     filter_name='TC_filter_H.dat', spectype="A0V")

point_src3 = simcado.source.stars(mags=magnitudes, x=x_pos, y=y_pos+14, 
                                     filter_name='TC_filter_H.dat', spectype="A0V")

point_src4 = simcado.source.stars(mags=magnitudes, x=x_pos, y=y_pos+21, 
                                     filter_name='TC_filter_H.dat', spectype="A0V")

src = point_src1 + point_src2 + point_src3 + point_src4


# Doing the simulation
t_exp = 100 # s of exposure time

sim_hdu = simcado.run(src, OBS_EXPTIME=t_exp, detector_layout="full", 
                  FPA_LINEARITY_CURVE=None, filter_name='TC_filter_H.dat',
                  SCOPE_PSF_FILE="AnisoCADO_SCAO_FVPSF_4mas_EsoMedian_20190328.fits")


# writing to disk
sim_hdu.writeto("random_stars.fits", overwrite=True)

