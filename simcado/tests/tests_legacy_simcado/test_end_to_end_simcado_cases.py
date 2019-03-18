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
    def see_what_happens(self):
        cmd = sim.UserCommands(sim_data_dir="C:/Work/Legacy_SimCADO_data/")


