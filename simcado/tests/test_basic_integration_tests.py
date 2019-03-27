import pytest
from astropy.io import fits
import simcado as sim


# class TestBasicUseCases:
#     def test_simple_run(self):
#         src = sim.source.cluster()
#         hdu = sim.run(src, OBS_EXPTIME=3600, SCOPE_PSF_FILE="PSF_MCAO.fits",
#                       filter_name="J", mode="zoom")
#         assert isinstance(hdu, fits.HDUList)
#
#


class Test