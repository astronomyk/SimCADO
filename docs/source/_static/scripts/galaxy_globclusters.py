#!/usr/bin/env python

"""
Generate a galaxy (exponential profile, r_eff=0.3", 50deg inclination),
with integrated magnitude K_AB=24. Put a few point sources (globular
cluster progenitors) around it with magnitudes in the range K_AB=26-28mag.
Set DIT=20sec and NDIT=6 (might need to increase this to see anything).
And take 4 exposures with the object shifted by ~1arcsec each time (e.g.
in a square). Use a MAORY-like PSF.

For simplicity here I'm doing one DIT with total integration time of 8000s ~ 2h
The file containing the MAORY-like PSF PSF_MCAO_Ks_Strehl40.fits must be placed
where SimCADO can find it, either with the full path or placing it in the SIM_DATA_DIR

see Notebook for more explanations and examples.
"""

import simcado
import matplotlib.pyplot as plt
import numpy as np
print(simcado.__data_dir__)

# Utility functions to calculate ellipticity from the inclination and vice versa
def inclination2ellipticity(inclination, q=0.2):
    """
    calculate galaxy ellipticity from the inclination (in degrees!)
    q is the true axis ratio with default q=0.2 as in Heidmann 1972
    """
    i = np.deg2rad(inclination)
    btoa = np.sqrt((1 - q**2) * np.cos(i)**2 + q**2)
    ell = 1 - btoa
    return ell

def ellipticity2inclination(ellipticity, q=0.2):
    """
    calculate galaxy inclination (in degrees) from the ellipticity measured on images
    q is the true axis ratio with default q=0.2 as in Heidmann 1972
    """
    btoa = 1 - ellipticity
    i =  np.arccos(np.sqrt( (btoa**2 - q**2)/(1 - q**2) ))
    inclination = np.rad2deg(i)
    return inclination


# Galaxy parameters
magAB = 22
magVega = magAB - 1.85
Reff = 0.3               # arcsec
n = 1                    # sersic index for a disk galaxy
inclination = 50         # deg
ellipticity = inclination2ellipticity(inclination)

# point source parameters
magnitudes = np.random.uniform(low=26, high=28, size=10) - 1.85
x_pos = np.random.normal(scale=Reff, size=10)
y_pos = np.random.normal(scale=Reff, size=10)

# Exposure times
t_exp = 400
NDIT  = 20
t_exp = t_exp*NDIT

# Doing four dithers in a square pattern
x_shifts = [1,  1, -1, -1]   # in arcsec
y_shifts = [1, -1,  1, -1]

n=0
for xs, ys in zip(x_shifts, y_shifts):
    print("creating sources at positions=", (xs, ys), "arcsec")
    gal_src = simcado.source.elliptical(Reff, 0.004, magnitude=magVega, n=1,
                                        x_offset= xs, y_offset=ys,
                                        filter_name='TC_filter_Ks.dat', spectrum="spiral",
                                        ellipticity=ellipticity)

    point_src = simcado.source.stars(mags=magnitudes, x=x_pos+xs, y=y_pos+ys,
                                     filter_name='TC_filter_Ks.dat', spectype="A0V")

    combined_src = gal_src + point_src

    hdu = simcado.run(combined_src, OBS_DIT=t_exp, detector_layout="small",
                  FPA_LINEARITY_CURVE=None, filter_name='TC_filter_Ks.dat',
                  SCOPE_PSF_FILE="PSF_MCAO_Ks_Strehl40.fits")

    filename = "globclustersKs_" + "dither" + str(n+1) + ".fits"
    hdu.writeto(filename, overwrite=True) # writing to fits
    n = n+1
