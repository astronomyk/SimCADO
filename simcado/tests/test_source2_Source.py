# actually for Source2
import pytest

import os
import inspect

import numpy as np

from astropy.io import fits
from astropy.io import ascii as ioascii
from astropy.table import Table
from astropy import units as u

from synphot import SourceSpectrum
from synphot.models import Empirical1D

import simcado as sim
from simcado.utils import convert_table_comments_to_dict
from simcado.source.source2 import Source


def mock_dir():
    cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
    rel_dirname = "mocks/sources/"
    return os.path.abspath(os.path.join(cur_dirname, rel_dirname))


MOCK_DIR = mock_dir()
sim.rc.__search_path__.insert(0, MOCK_DIR)


@pytest.fixture(scope="module")
def input_files():
    filenames = ["test_image.fits", "test_table.fits", "test_table.tbl",
                 "test_spectrum_Flam.dat", "test_spectrum_photlam.dat"]
    filenames = [os.path.join(MOCK_DIR, fname) for fname in filenames]
    return filenames


@pytest.fixture(scope="module")
def input_imagehdu():
    filenames = ["test_image.fits"]
    filenames = [os.path.join(MOCK_DIR, fname) for fname in filenames]
    hdu_handle = fits.open(filenames[0])

    return hdu_handle


@pytest.fixture(scope="module")
def input_tables():
    filenames = ["test_table.fits", "test_table.tbl"]
    filenames = [os.path.join(MOCK_DIR, fname) for fname in filenames]
    tbls = []
    tbls += [Table.read(filenames[0])]
    tbls += [Table.read(filenames[1], format="ascii.basic")]
    tbls[1].meta.update(convert_table_comments_to_dict(tbls[1]))

    return tbls


@pytest.fixture(scope="module")
def input_spectra():
    filenames = ["test_spectrum_photlam.dat", "test_spectrum_Flam.dat"]
    filenames = [os.path.join(MOCK_DIR, fname) for fname in filenames]
    tbls = [ioascii.read(fname) for fname in filenames]
    specs = []
    for tbl in tbls:
        tbl.meta.update(convert_table_comments_to_dict(tbl))
        wave = tbl["wavelength"] * u.Unit(tbl.meta["wavelength_unit"])
        flux = tbl["flux"] * u.Unit(tbl.meta["flux_unit"])
        specs += [SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)]

    return specs


@pytest.mark.usefixtures("input_files", "input_imagehdu", "input_tables",
                         "input_spectra")
class TestSourceInit:
    def test_initialises_with_nothing(self):
        src = Source()
        assert isinstance(src, Source)

    @pytest.mark.parametrize("ii", [0, 1, 2])
    def test_initialises_with_table_and_2_spectrum(self, ii,
                                                   input_tables,
                                                   input_spectra):
        table = input_tables[0]
        src = Source(table=table, spectra=input_spectra)
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.positions[0], Table)

    # @pytest.mark.parametrize("ii", [0, 1, 2])
    # def test_initialises_with_filename_and_spectrum(self, ii, input_files,
    #                                                 input_spectra):
    #     fname = input_files[ii]
    #     src = Source(filename=fname, spectra=input_spectra)
    #     assert isinstance(src, Source)
    #     assert isinstance(src.spectra[0], SourceSpectrum)
    #     assert isinstance(src.positions[0], Table)
    #
