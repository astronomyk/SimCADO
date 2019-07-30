#!/usr/bin/env python

"""
Use an image of an extended source generated externally. Import it, and 
scale it to a given integrated K-band magnitude. Observe it with a SCAO 
PSF for a 16mag on-axis guide star.

Using Halpha image of HCG04 from Eigenthaler et al. 2015 
http://adsabs.harvard.edu/abs/2015MNRAS.451.2793E, because it has fewer contaminating nearby objects. 

Main parameters Distance = 106.36 Mpc Scale ~ 515 pc/arcsec pixel_scale=0.145" -> 74.6 pc/pixel. MICADO will provide ~33 pc/pixel at z=2. 
Image will be undersampled, but for a factor ~2. 

Data for this example can be found here 
https://github.com/astronomyk/SimCADO/tree/master/docs/source/_static/images/HCG04_Halpha_fluxcalibrated.fits
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import simcado
print(simcado.__data_dir__)

original = fits.getdata("HCG04_Halpha_fluxcalibrated.fits")
galaxy_data=original[1100:1800,1100:1800] # cutting the edges


z = 2
magnitude = 20 # Vega system, probably the magnitude of a L* galaxy at z~2
               # Pozzetti, 2003, Cirasuolo 2007, Gobat 2011

ec = simcado.source.redshift_SED(z, "spiral", mag=20, filter_name='TC_filter_Ks.dat')
lam, spec = ec.lam, ec.val

galaxy_src = simcado.source.source_from_image(galaxy_data, lam, spec, 
    plate_scale=0.004, flux_threshold=0.01)

# plate scale change the apparent size of the galaxy, here is 1:1
# flux_threshold rejects all pixels under that value

t_exp = 3600*2  # 2h exp time

sim_galaxy = simcado.run(galaxy_src, OBS_EXPTIME=t_exp, 
                         INST_FILTER_TC="TC_filter_Ks.dat", SCOPE_PSF_FILE="PSF_SCAO.fits",
                         FPA_LINEARITY_CURVE=None)

sim_galaxy.writeto("sim_distant_galaxy.fits", overwrite=True)