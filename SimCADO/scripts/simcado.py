#!/usr/bin/env python

##################
# SimCADO script #
##################
'''SimCADO script'''

########   Imports    ##########
import os
import numpy as np
from astropy import units as u
from astropy.io import fits
from sys import argv

import SimCADO as sm

def run_simcado():
    """Function to run simcado"""

    ### Config ###
    config_dict = sm.read_config("simcado.default")
    if len(argv) > 1:
        user_dict = sm.read_config(argv[1])
        config_dict.update(user_dict)

    # Alternativ:
    # config_dict = sm.update_config(argv[1], config_dict)

    ### Pull in TC files ###
    if config_dict["ATMO_USE_SKYCALC_TC"] == 'yes' \
       and os.path.exists(config_dict["ATMO_SKYCALC_FILE"]):
        tc_atmo = sm.Throughput.from_skycalc(config_dict["ATMO_SKYCALC_FILE"])
    else:
        tc_atmo = sm.Throughput.from_ascii(config_dict["ATMO_TC"])

    tc_filter = sm.Throughput.from_ascii(config_dict["INST_TC_FILTER"],
                                         Type="Filter")
    tc_mirror = sm.Throughput.from_ascii(config_dict["INST_TC_MIRROR"],
                                         Type="Mirror")
    tc_detector = sm.Throughput.from_ascii(config_dict["FPA_QE"])

    n_mirror = int(config_dict["SCOPE_NUM_MIRRORS"]) + \
               int(config_dict["INST_NUM_MIRRORS"])

    surfaces = []

    for i in [tc_atmo] + [tc_mirror] * n_mirror + [tc_filter, tc_detector]:
        if i is not None:
            surfaces.append(i)

    tc_system = sm.Throughput.from_throughput(surfaces)

    ## TODO
    ## Next, Kieran determines a wavelength grid, defined by values
    ## obs_lam_min, obs_lam_max
    ## obs_lam_res
    ## obs_pix_res
    obs_pix_res  = config_dict["OBS_PIXEL_SCALE"] * u.mas
    ## lam_bin_edges
    ## lam_bin_centre
    

    ### Generate the PSFs ###

    ao_effect = (100. - float(config_dict["SCOPE_AO_EFFECTIVENESS"])) * 6 * u.mas

    psf_telescope = sm.PSFCube.gen_cube(lam_bin_centre, res=obs_pix_res,
                                        kernel='airy')

    psf_atmo = sm.PSFCube.gen_cube(lam_bin_centre, res=obs_pix_res,
                                   kernel='gauss', fwhm=ao_effect)

    psf_jitter = sm.PSFCube.gen_cube(lam_bin_centre, res=obs_pix_res,
                                     kernel='gauss',
                                     fwhm=float(config_dict["SCOPE_JITTER_FWHM"]))
    
    psf_adc = sm.PSF_Cube.gen_adc(lam_bin_centre, config_dict, res=obs_pix_res)

    ## Kieran's system_psf works differently (explicit use of collapse())
    psf_system = sm.PSFCube.gen_from_list([psf_telescope, psf_adc,
                                           psf_atmo, psf_jitter])

    ## Get source cube





