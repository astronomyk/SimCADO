#!/usr/bin/env python

import simcado as sim
sim.get_extras()
mcao_fv_psf_url = "http://www.univie.ac.at/simcado/InstPkgSvr/psfs/MAORY_MCAO_FVPSF_4mas_20181203.fits"
sim.utils.download_file(mcao_fv_psf_url, sim.__data_dir__)
