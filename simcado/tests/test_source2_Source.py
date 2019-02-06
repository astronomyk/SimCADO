# actually for Source2
import pytest

import os
import inspect

import numpy as np

from astropy.io import fits
from astropy.io import ascii as ioascii
from astropy.table import Table
from astropy import units as u
from astropy import wcs

from synphot import SourceSpectrum, SpectralElement, Observation
from synphot.models import Empirical1D
from synphot.units import PHOTLAM

import simcado as sim
from simcado.utils import convert_table_comments_to_dict
from simcado.source.source2 import Source
from simcado.source import source2 as src2


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
def input_hdulist():
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


@pytest.mark.usefixtures("input_files", "input_hdulist", "input_tables",
                         "input_spectra")
class TestSourceInit:
    def test_initialises_with_nothing(self):
        src = Source()
        assert isinstance(src, Source)

    @pytest.mark.parametrize("ii", [0, 1])
    def test_initialises_with_table_and_2_spectrum(self, ii,
                                                   input_tables,
                                                   input_spectra):
        table = input_tables[ii]
        src = Source(table=table, spectra=input_spectra)
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.positions[0], Table)

    def test_initialises_with_image_and_1_spectum(self, input_hdulist,
                                                  input_spectra):
        src = Source(image=input_hdulist[0], spectra=input_spectra)
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.positions[0], fits.PrimaryHDU)

    @pytest.mark.parametrize("ii, dtype",
                             [(0, fits.ImageHDU),
                              (1, Table),
                              (2, Table)])
    def test_initialises_with_filename_and_spectrum(self, ii, dtype,
                                                    input_files, input_spectra):
        fname = input_files[ii]
        src = Source(filename=fname, spectra=input_spectra)
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.positions[0], dtype)

    def test_initialised_with_old_style_arrays(self):
        x, y = [0, 1], [0, -1]
        ref, weight = [0, 0], [1, 10]
        lam = np.linspace(0.5, 2.5, 11) * u.um
        spectra = np.ones(11) * PHOTLAM
        src = Source(x=x, y=y, ref=ref, weight=weight, lam=lam, spectra=spectra)
        assert isinstance(src, Source)
        assert isinstance(src.spectra[0], SourceSpectrum)
        assert isinstance(src.positions[0], Table)


class TestSourceAddition:
    def test_ref_column_always_references_correct_spectrum(self):
        pass


@pytest.mark.usefixtures("input_spectra")
class TestPhotonsInRange:
    @pytest.mark.parametrize("ii, n_ph",
                             [(0, 50),
                              (1, 8e12)])
    def test_returns_correct_number_of_photons_for_one_spectrum(self, ii, n_ph,
                                                                input_spectra):
        spec = input_spectra[ii]
        counts = src2.photons_in_range([spec], 1, 2, nbins=100)
        assert np.isclose(counts.value, n_ph, rtol=2e-3)

    def test_returns_ones_for_unity_spectrum(self):
        flux = np.ones(11) * u.Unit("ph s-1 m-2 um-1")
        wave = np.linspace(1, 2, 11) * u.um
        spec = SourceSpectrum(Empirical1D, points=wave, lookup_table=flux)
        counts = src2.photons_in_range([spec], 1*u.um, 2*u.um)
        assert np.isclose(counts.value, 1, rtol=2e-4)


class TestMakeImageFromTable:
    def test_returned_object_is_image_hdu(self):
        hdu = src2.make_imagehdu_from_table(x=[0], y=[0], flux=[1])
        assert isinstance(hdu, fits.ImageHDU)

    def test_imagehdu_has_wcs(self):
        hdu = src2.make_imagehdu_from_table(x=[0], y=[0], flux=[1])
        wcs_keys = ["CRPIX1", "CRVAL1", "CDELT1", "CUNIT1"]
        assert np.all([key in hdu.header.keys() for key in wcs_keys])

    def test_stars_are_at_correct_position_for_no_subpixel_accuracy(self):
        # weird behaviour for exactly integer coords e.g. (0, 0)
        x = np.linspace(-1.0001, 1.0001, 9)*u.arcsec
        y = np.linspace(-1.0001, 1.0001, 9)*u.arcsec
        flux = np.ones(len(x))
        hdu = src2.make_imagehdu_from_table(x=x, y=y, flux=flux,
                                           pix_scale=0.25*u.arcsec)
        the_wcs = wcs.WCS(hdu)
        yy, xx = the_wcs.wcs_world2pix(y.to(u.deg), x.to(u.deg), 1)
        xx = np.floor(xx).astype(int)
        yy = np.floor(yy).astype(int)
        assert np.all(hdu.data[xx, yy])

    @pytest.mark.parametrize("ii", np.arange(20))
    def test_wcs_returns_correct_pixel_values_for_random_coords(self, ii):
        x = np.random.random(11)*u.arcsec
        y = np.random.random(11)*u.arcsec
        flux = np.ones(len(x))
        hdu = src2.make_imagehdu_from_table(x=x, y=y, flux=flux,
                                           pix_scale=0.1*u.arcsec)
        the_wcs = wcs.WCS(hdu)
        yy, xx = the_wcs.wcs_world2pix(y.to(u.deg), x.to(u.deg), 1)
        xx = xx.astype(int)
        yy = yy.astype(int)
        assert np.all(hdu.data[xx, yy])

        # When plotting the image vs the scatter plot
        # The dots match up with a 0.5 px shift, When we intruduce the shift to
        # the test, the dots are on top of the image, but the test fails
        # when using
        # the_wcs.wcs.crpix = [0.5, 0.5]
        # Returning to normal the test passes when (albeit with an image offset)
        # the_wcs.wcs.crpix = [0., 0.]
        #
        # print(hdu.data, flush=True)
        # print(x, y, flush=True)
        # import matplotlib.pyplot as plt
        # plt.subplot(projection=the_wcs)
        # plt.scatter(y*10, x*10)
        # plt.scatter(yy , xx)
        # plt.imshow(hdu.data)
        # plt.show()


class TestSourceImageInRange:
    def test_that_it_does_what_it_should(self):
        pass
