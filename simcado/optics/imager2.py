import os
import warnings
from copy import deepcopy

from astropy.io import ascii as ioascii
from astropy.table import Table, vstack
from astropy import units as u

from synphot import SpectralElement
from synphot.models import Empirical1D


class Imager:
    def __init__(self, cmds=None):
        self.cmds = deepcopy(cmds)

    @property
    def surfaces(self):
        files = [self.cmds["SCOPE_MIRROR_LIST"],
                 self.cmds["INST_MIRROR_AO_LIST"],
                 self.cmds["INST_MIRROR_LIST"]]
        surf_tbl = make_surfaces_table(filenames=files)

        return surf_tbl



def make_surfaces_table(filenames=()):

    if isinstance(filenames, str):
        filenames = [filenames]

    if len(filenames) == 0:
        warnings.warn("'filenames' was empty. No tables were read")

    lst = [ioascii.read(fname) for fname in filenames
           if fname is not None and os.path.exists(fname)]

    if len(lst) > 0:
        tbl = vstack(lst)
    else:
        tbl = Table()
    return tbl


def make_spectral_curve_dict(filenames=()):
    pass


def import_spectral_curve_from_file(filename,
                                    wave_name="wavelength", wave_unit="um",
                                    val_name="transmission", val_unit=""):

    if not os.path.exists(filename):
        raise ValueError("{} doesn't exist".format(filename))

    tbl = ioascii.read(filename)

    for col_name in [wave_name, val_name]:
        if col_name not in tbl.colnames:
            raise ValueError("{} column does not exist in {}".format(col_name,
                                                                     filename))
    spec_curve = SpectralElement(Empirical1D,
                                 points=tbl[wave_name]*u.Unit(wave_unit),
                                 lookup_table=tbl[val_name]*u.Unit(val_unit))
    return spec_curve


