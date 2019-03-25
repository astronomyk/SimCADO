import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

import simcado as sim
import simcado.rc

PLOTS = False
PKG_MICADO = "C:/Work/irdb/_Legacy_packages/MICADO/"
FVPSF_PATH = "C:/Work/Legacy_SimCADO_data/MAORY_SCAO_FVPSF_4mas_20181203.fits"


class TestSimDataDir:
    def test_user_commands_defaults_to_installation_directory(self):
        assert len(simcado.rc.__search_path__) > 1


class TestNormalSimcadoUse:
    def test_basic_cluster_example(self):
        src = sim.source.cluster(mass=1e4, distance=50e3)
        hdu = sim.run(src, sim_data_dir=PKG_MICADO)
        if PLOTS:
            plt.imshow(hdu[0].data, norm=LogNorm())
            plt.show()


class TestFVPSFsWithSimCADO:
    def test_see_what_happens(self):
        cmd = sim.UserCommands(sim_data_dir=PKG_MICADO)
        cmd["SCOPE_PSF_FILE"] = FVPSF_PATH
        opt = sim.OpticalTrain(cmd)
        assert isinstance(opt.psf, sim.psf.FieldVaryingPSF)

    def test_reading_in_fv_psf_with_legacy_functions(self):
        cmd = sim.UserCommands(sim_data_dir=PKG_MICADO)
        cmd["SCOPE_PSF_FILE"] = FVPSF_PATH

        from simcado.psf import FieldVaryingPSF
        psf = FieldVaryingPSF(filename=cmd["SCOPE_PSF_FILE"])
        print(psf.info)

        if PLOTS:
            plt.imshow(psf.strehl_imagehdu.data.T, origin="l")
            plt.show()

    def what_happens_when_passed_to_source(self):
        cmd = sim.UserCommands(sim_data_dir="C:/Work/Legacy_SimCADO_data/")
        cmd["FPA_LINEARITY_CURVE"] = None
        cmd["SIM_USE_FILTER_LAM"] = "no"
        cmd["SIM_LAM_MIN"] = 1.9
        cmd["SIM_LAM_MAX"] = 2.4
        cmd["OBS_EXPTIME"] = 3600
        cmd["SCOPE_PSF_FILE"] = "MAORY_SCAO_FVPSF_4mas_20181203.fits"
        cmd["FPA_CHIP_LAYOUT"] = "full"

        opt = sim.OpticalTrain(cmd)
        fpa = sim.Detector(cmd, small_fov=False)

        src = sim.source.cluster(mass=1e4, distance=20e3)

        src.x = -24 + 48.*np.random.random(len(src.x))
        src.y = -24 + 48.*np.random.random(len(src.y))
        src.apply_optical_train(opt, fpa)
        hdu = fpa.read_out()
        hdu.writeto("E:/fv_psf.fits", clobber=True)

        if PLOTS:
            plt.imshow(hdu[0].data.T, origin="l")
            plt.show()


class TestPoorMansFOV:
    def initialised_with_a_chip(self):
        from simcado.fv_psf import PoorMansFOV

        cmd = sim.UserCommands(sim_data_dir=PKG_MICADO)
        cmd["FPA_CHIP_LAYOUT"] = "center"
        fpa = sim.Detector(cmd, small_fov=False)
        fov = PoorMansFOV(fpa.chips[0], cmd.lam_bin_edges[0],
                          cmd.lam_bin_edges[-1])
        assert isinstance(fov, PoorMansFOV)


class TestCaseStudiesForFVPSFs:
    def grid_of_stars(self):
        import time
        start = time.time()

        cmd = sim.UserCommands(sim_data_dir="C:/Work/Legacy_SimCADO_data/")
        cmd["FPA_LINEARITY_CURVE"] = None
        cmd["SIM_USE_FILTER_LAM"] = "no"
        cmd["INST_FILTER_TC"] = "TC_filter_J.dat"
        cmd["SIM_LAM_MIN"] = 1.1
        cmd["SIM_LAM_MAX"] = 1.35
        cmd["OBS_EXPTIME"] = 3600
        mcao = "C:/Work/irdb/_PSFs/MAORY_MCAO_FVPSF_4mas_20181203.fits"
        scao = "C:/Work/irdb/_PSFs/AnisoCADO_SCAO_FVPSF_4mas_1024_20190321.fits"
        cmd["SCOPE_PSF_FILE"] = scao
        cmd["FPA_CHIP_LAYOUT"] = "centre"
        cmd["OBS_SCAO_NGS_OFFSET_X"] = 0
        cmd["OBS_SCAO_NGS_OFFSET_Y"] = 0

        opt = sim.OpticalTrain(cmd)
        fpa = sim.Detector(cmd, small_fov=False)

        # src = sim.source.star_grid(900, 15, 15.1, separation=2)
        src = sim.source.cluster(mass=1E4, distance=20000, half_light_radius=0.5)

        src.apply_optical_train(opt, fpa)
        hdu = fpa.read_out()
        hdu.writeto("E:/test_psf_aniso_J_centre.fits", clobber=True)

        end = time.time()
        print("Time elapsed: {} sec".format(end - start))

        # import matplotlib.pyplot as plt
        # from matplotlib.colors import LogNorm
        # plt.imshow(hdu[0].data.T, norm=LogNorm())
        # plt.show()
