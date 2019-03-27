import os
import pytest
from astropy.io import fits
import simcado as sim


# class TestBasicUseCases:
#     def test_simple_run(self):
#         src = sim.source.cluster()
#         hdu = sim.run(src, OBS_EXPTIME=3600, SCOPE_PSF_FILE="PSF_MCAO.fits",
#                       filter_name="J", mode="zoom")
#         assert isinstance(hdu, fits.HDUList)


class TestTravisDownloadsGetExtras:
    def test_find_psf_file_with_find_file(self):
        dirpath = "C:\Work\irdb\_Legacy_packages\MICADO"
        filename = "PSF_SCAO.fits"
        filepath = sim.utils.find_file(filename, path=[dirpath])
        print(filepath, os.path.join(dirpath, filename))

        assert os.path.exists(filepath) is not None
