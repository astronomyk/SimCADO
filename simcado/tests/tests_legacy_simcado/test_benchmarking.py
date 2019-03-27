# import time
# import numpy as np
#
# from matplotlib import pyplot as plt
# from matplotlib.colors import LogNorm
#
# import simcado as sim
#
# PLOTS = False
#
# # To run the test, add "test_" back into the names of the functions
#
#
# class TestConstantPSFs:
#     def basic_cluster_example_1024_px_window(self):
#         start = time.time()
#
#         src = sim.source.cluster(mass=1e4, distance=20e3, half_light_radius=0.1)
#         hdu = sim.run(src, sim_data_dir="C:/Work/Legacy_SimCADO_data/")
#
#         hdu.writeto("E:/test_psf.fits", clobber=True)
#
#         if PLOTS:
#             plt.imshow(hdu[0].data.T, origin="l")
#             plt.show()
#
#         end = time.time()
#         print("Time elapsed: {} sec".format(end - start))
#
#     def basic_cluster_example_central_detector(self):
#         start = time.time()
#
#         src = sim.source.cluster(mass=1e4, distance=20e3)
#         hdu = sim.run(src, sim_data_dir="C:/Work/Legacy_SimCADO_data/",
#                       detector_layout="centre")
#
#         hdu.writeto("E:/test_psf.fits", clobber=True)
#
#         end = time.time()
#         print("Time elapsed: {} sec".format(end - start))
#
#     def basic_cluster_example_full_fpa(self):
#         start = time.time()
#
#         src = sim.source.cluster(mass=1e4, distance=20e3)
#         hdu = sim.run(src, sim_data_dir="C:/Work/Legacy_SimCADO_data/",
#                       detector_layout="full")
#
#         hdu.writeto("E:/test_psf.fits", clobber=True)
#
#         end = time.time()
#         print("Time elapsed: {} sec".format(end - start))
#
#
# class TestCaseFVPSFs:
#
#     def basic_cluster_example_1024_px_window(self):
#         start = time.time()
#
#         cmd = sim.UserCommands(
#             sim_data_dir="C:/Work/Legacy_SimCADO_data/")
#         cmd["FPA_LINEARITY_CURVE"] = None
#         cmd["SIM_USE_FILTER_LAM"] = "no"
#         cmd["INST_FILTER_TC"] = "TC_filter_J.dat"
#         cmd["SIM_LAM_MIN"] = 1.0
#         cmd["SIM_LAM_MAX"] = 1.25
#         cmd["OBS_EXPTIME"] = 3600
#         cmd["SCOPE_PSF_FILE"] = "MAORY_SCAO_FVPSF_4mas_20181203.fits"
#         cmd["FPA_CHIP_LAYOUT"] = "small"
#
#         src = sim.source.cluster(mass=1e4, distance=50e3)
#         hdu = sim.run(src, cmds=cmd,
#                       sim_data_dir="C:/Work/Legacy_SimCADO_data/")
#
#         hdu.writeto("E:/test_psf.fits", clobber=True)
#
#         end = time.time()
#         print("Time elapsed: {} sec".format(end - start))
#
#     def grid_of_stars_full_fpa(self):
#         start = time.time()
#
#         cmd = sim.UserCommands(
#             sim_data_dir="C:/Work/Legacy_SimCADO_data/")
#         cmd["FPA_LINEARITY_CURVE"] = None
#         cmd["SIM_USE_FILTER_LAM"] = "no"
#         cmd["INST_FILTER_TC"] = "TC_filter_J.dat"
#         cmd["SIM_LAM_MIN"] = 1.0
#         cmd["SIM_LAM_MAX"] = 1.25
#         cmd["OBS_EXPTIME"] = 3600
#         cmd["SCOPE_PSF_FILE"] = "MAORY_SCAO_FVPSF_4mas_20181203.fits"
#         cmd["FPA_CHIP_LAYOUT"] = "full"
#
#         opt = sim.OpticalTrain(cmd)
#         fpa = sim.Detector(cmd, small_fov=False)
#
#         src = sim.source.star_grid(900, 15, 15.1, separation=2)
#         src.apply_optical_train(opt, fpa)
#         hdu = fpa.read_out()
#
#         hdu.writeto("E:/test_psf.fits", clobber=True)
#
#         end = time.time()
#         print("Time elapsed: {} sec".format(end - start))
#
#
