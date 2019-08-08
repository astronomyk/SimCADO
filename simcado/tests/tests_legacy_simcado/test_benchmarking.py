import os
import time

import numpy as np
from astropy.io import fits

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

import simcado as sim

PLOTS = False

if "USERNAME" in os.environ and os.environ["USERNAME"] == "Kieran":
    sim.__search_path__.insert(0, "C:/Work/irdb/_Legacy_packages/MICADO")
    sim.__search_path__.insert(1, "C:/Work/Legacy_SimCADO_data/")


class TestConstantPSFs:
    def test_basic_cluster_example_1024_px_window(self):
        start = time.time()

        src = sim.source.cluster(mass=1e4, distance=20e3, half_light_radius=0.1)
        hdu = sim.run(src)
        # hdu.writeto("E:/test_psf.fits", clobber=True)

        if PLOTS:
            plt.imshow(hdu[0].data.T, origin="l")
            plt.show()

        end = time.time()
        print("Time elapsed: {} sec".format(end - start))

        assert isinstance(hdu, fits.HDUList)

    def test_basic_cluster_example_central_detector(self):
        start = time.time()

        src = sim.source.cluster(mass=1e4, distance=20e3)
        hdu = sim.run(src, detector_layout="centre")
        # hdu.writeto("E:/test_psf.fits", clobber=True)

        end = time.time()
        print("Time elapsed: {} sec".format(end - start))

        assert isinstance(hdu, fits.HDUList)

    def test_basic_cluster_example_full_fpa(self):
        start = time.time()

        src = sim.source.cluster(mass=1e4, distance=20e3)
        hdu = sim.run(src, detector_layout="full")

        # hdu.writeto("E:/test_psf.fits", clobber=True)

        end = time.time()
        print("Time elapsed: {} sec".format(end - start))

        assert isinstance(hdu, fits.HDUList)


class TestCaseFVPSFs:
    def test_basic_cluster_example_1024_px_window(self):
        start = time.time()

        cmd = sim.UserCommands()
        cmd["FPA_LINEARITY_CURVE"] = None
        cmd["SIM_USE_FILTER_LAM"] = "no"
        cmd["INST_FILTER_TC"] = "TC_filter_J.dat"
        cmd["SIM_LAM_MIN"] = 1.0
        cmd["SIM_LAM_MAX"] = 1.25
        cmd["OBS_DIT"] = 3600
        cmd["SCOPE_PSF_FILE"] = "MAORY_MCAO_FVPSF_4mas_20181203.fits"
        cmd["FPA_CHIP_LAYOUT"] = "small"

        src = sim.source.cluster(mass=1e4, distance=50e3)
        hdu = sim.run(src, cmds=cmd)
        # hdu.writeto("E:/test_psf.fits", clobber=True)

        end = time.time()
        print("Time elapsed: {} sec".format(end - start))

    def grid_of_stars_full_fpa(self):
        start = time.time()

        cmd = sim.UserCommands()
        cmd["FPA_LINEARITY_CURVE"] = None
        cmd["SIM_USE_FILTER_LAM"] = "no"
        cmd["INST_FILTER_TC"] = "TC_filter_J.dat"
        cmd["SIM_LAM_MIN"] = 1.0
        cmd["SIM_LAM_MAX"] = 1.25
        cmd["OBS_EXPTIME"] = 3600
        cmd["SCOPE_PSF_FILE"] = "MAORY_SCAO_FVPSF_4mas_20181203.fits"
        cmd["FPA_CHIP_LAYOUT"] = "full"

        opt = sim.OpticalTrain(cmd)
        fpa = sim.Detector(cmd, small_fov=False)

        src = sim.source.star_grid(900, 15, 15.1, separation=2)
        src.apply_optical_train(opt, fpa)
        hdu = fpa.read_out()
        #hdu.writeto("E:/test_psf.fits", clobber=True)

        end = time.time()
        print("Time elapsed: {} sec".format(end - start))
