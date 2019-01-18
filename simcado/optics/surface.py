import os

from astropy import table
import synphot

class SpectralSurface:
    def __init__(self, filename=None, **kwargs):
        self.params = {""}
        self.params.update(kwargs)

        pass

    @property
    def transmission(self):
        pass

    @property
    def emission(self):
        pass

    @property
    def transmission(self):
        pass



# from astropy.io import ascii
# import glob
# import os
#
#
# if __name__ == "__main__":
#
#     dirname = "E:\Work\SimCADO_inst_pkgs\MICADO"
#     for fname in glob.glob(os.path.join(dirname, "TC_filter*.*")):
#
#         ascii.read(fname)
#         print(fname)
