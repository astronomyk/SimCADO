import numpy as np

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

import simcado as sim

PLOTS = True


class TestSimDataDir:
    def test_user_commands_defaults_to_installation_directory(self):
        print(sim.__data_dir__)


class TestNormalSimcadoUse:
    def test_basic_cluster_example(self):
        src = sim.source.cluster(mass=1e4, distance=50e3)
        hdu = sim.run(src, sim_data_dir="C:/Work/Legacy_SimCADO_data/")
        if PLOTS:
            plt.imshow(hdu[0].data, norm=LogNorm())
            plt.show()


class TestFVPSFsWithSimCADO:
    def test_see_what_happens(self):
        cmd = sim.UserCommands(sim_data_dir="C:/Work/Legacy_SimCADO_data/")
        cmd["SCOPE_PSF_FILE"] = "MAORY_SCAO_FVPSF_4mas_20181203.fits"
        opt = sim.OpticalTrain(cmd)
        print(opt.psf)

    def test_reading_in_fv_psf_with_legacy_functions(self):
        cmd = sim.UserCommands(sim_data_dir="C:/Work/Legacy_SimCADO_data/")
        cmd["SCOPE_PSF_FILE"] = "MAORY_SCAO_FVPSF_4mas_20181203.fits"
        from simcado.psf import FieldVaryingPSF
        psf = FieldVaryingPSF(filename=cmd["SCOPE_PSF_FILE"])
        print(psf.info)

        if PLOTS:
            plt.imshow(psf.strehl_imagehdu.data.T, origin="l")
            plt.show()

    def test_what_happens_when_passed_to_source(self):
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
        hdu.writeto("E:/test_fv_psf.fits", clobber=True)

        plt.imshow(hdu[0].data.T, origin="l")
        plt.show()


class TestPoorMansFOV:
    def test_initialised_with_a_chip(self):
        from simcado.fv_psf import PoorMansFOV

        cmd = sim.UserCommands(sim_data_dir="C:/Work/Legacy_SimCADO_data/")
        cmd["FPA_CHIP_LAYOUT"] = "center"
        fpa = sim.Detector(cmd, small_fov=False)
        fov = PoorMansFOV(fpa.chips[0], cmd.lam_bin_edges[0],
                          cmd.lam_bin_edges[-1])
        assert isinstance(fov, PoorMansFOV)

        print(dict(fov.hdu.header))
