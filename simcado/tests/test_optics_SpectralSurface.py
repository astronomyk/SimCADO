# 1 read in a table
# 2 compliment the table based on columns in file
# 3 have @property methods for: transmission, ermission, reflection

import pytest
import inspect
import os
import sys
import warnings

import numpy as np
from astropy.table import Table, Column
from astropy import units as u
from synphot import SpectralElement, SourceSpectrum

from simcado.optics import surface as opt_surf


@pytest.fixture(scope="module")
def input_tables():
    cur_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))
    rel_dirname = "../tests/mocks/MICADO_SCAO_WIDE/"
    abs_dirname = os.path.abspath(os.path.join(cur_dirname, rel_dirname))

    filenames = ["TER_dichroic.dat", "TC_filter_Ks.dat"]
    abs_paths = [os.path.join(abs_dirname, fname) for fname in filenames]

    return abs_paths


@pytest.fixture(scope="module")
def unity_flux():
    flux = np.ones(100)
    wave = np.logspace(-1, 1, 100) * u.um

    return flux, wave


@pytest.mark.usefixtures("input_tables")
class TestInit:
    def test_can_exist_with_no_input(self):
        srf = opt_surf.SpectralSurface()
        assert isinstance(srf, opt_surf.SpectralSurface)

    def test_reads_in_table_which_exists(self, input_tables):
        srf = opt_surf.SpectralSurface(filename=input_tables[0])
        assert isinstance(srf.table, Table)

    def test_returns_empty_table_if_path_is_bogus(self):
        srf = opt_surf.SpectralSurface(filename="bogus.txt")
        assert isinstance(srf, opt_surf.SpectralSurface)
        assert len(srf.table) == 0


@pytest.mark.usefixtures("input_tables")
class TestWavelengthProperty:
    def test_returns_quantity_array_from_file(self, input_tables):
        srf = opt_surf.SpectralSurface(filename=input_tables[0])
        assert isinstance(srf.wavelength, u.Quantity)
        assert srf.wavelength.unit == u.um

    def test_returns_quantity_array_from_file_with_no_unit(self, input_tables):
        srf = opt_surf.SpectralSurface(filename=input_tables[1])
        assert isinstance(srf.wavelength, u.Quantity)
        assert srf.wavelength.unit == u.um

    def test_returns_quantity_if_wavelength_overridden_by_list(self):
        srf = opt_surf.SpectralSurface(wavelength=[0.3, 3.0])
        assert np.all(srf.wavelength == [0.3, 3.0]*u.um)

    def test_returns_quantity_if_wavelength_unit_overridden_by_quantity(self):
        srf = opt_surf.SpectralSurface(wavelength=[0.3, 3.0],
                                       wavelength_unit=u.Angstrom)
        assert srf.meta["wavelength_unit"] == u.Angstrom
        assert np.all(srf.wavelength == [0.3, 3.0]*u.Angstrom)

    def test_returns_quantity_if_wavelength_unit_overridden_by_string(self):
        srf = opt_surf.SpectralSurface(wavelength=[0.3, 3.0],
                                       wavelength_unit="Angstrom")
        assert srf.meta["wavelength_unit"] == u.Angstrom
        assert np.all(srf.wavelength == [0.3, 3.0]*u.Angstrom)


@pytest.mark.usefixtures("input_tables")
class TestTransmissionProperty:
    def test_returns_synphot_object_array_from_file(self, input_tables):
        srf = opt_surf.SpectralSurface(filename=input_tables[0])
        assert isinstance(srf.transmission, SpectralElement)

    def test_returns_synphot_object_when_given_arrays(self):
        srf = opt_surf.SpectralSurface(wavelength=[0.3, 3.0]*u.um,
                                       transmission=[1, 1])
        assert isinstance(srf.transmission, SpectralElement)

    def test_returns_synphot_object_if_only_reflection_is_given(self):
        srf = opt_surf.SpectralSurface(wavelength=[0.3, 3.0]*u.um,
                                       reflection=[1, 1])
        assert isinstance(srf.transmission, SpectralElement)

    def test_returns_none_when_none_of_the_three_TER_info_is_given(self):
        srf = opt_surf.SpectralSurface(wavelength=[0.3, 3.0]*u.um)
        assert srf.transmission is None

    def test_returns_none_when_no_wavelength_info_given(self):
        srf = opt_surf.SpectralSurface(transmission=[1, 1])
        assert srf.transmission is None


class TestComplimentArray:
    @pytest.mark.parametrize("colname1, colname2, col1, col2, expected",
                             [("A", "B", [0.8]*u.um, [0.1]*u.um, [0.1]*u.um),
                              ("A", "B", [0.8]*u.um, None,       [0.2]*u.um),
                              ("A", "B", None,       [0.8]*u.um, [0.2]*u.um)])
    def test_the_right_answers_for_valid_input(self, colname1, colname2,
                                               col1, col2, expected):
        srf = opt_surf.SpectralSurface()
        srf.meta[colname1] = col1
        srf.meta[colname2] = col2
        col3 = srf._compliment_array(colname1, colname2)
        if sys.version_info.major >= 3:
            assert np.all(np.isclose(col3.data, expected.data))
            assert col3.unit == expected.unit
        else:
            warnings.warn("Data equality isn't tested for 2.7")
            assert col3.unit == expected.unit

    @pytest.mark.parametrize("colname1, colname2, col1, col2, expected",
                             [("A",     "B",      None, None, None)])
    def test_returns_none_for_none_input(self, colname1, colname2,
                                         col1, col2, expected):
        srf = opt_surf.SpectralSurface()
        srf.meta[colname1] = col1
        srf.meta[colname2] = col2
        col3 = srf._compliment_array(colname1, colname2)
        assert col3 is None

    @pytest.mark.parametrize("col2_arr, expected",
                             [([0.1],      [0.1]),
                              ([0.1, 0.1], [0.1, 0.1]),
                              ([0.1, 0.1, 0.1], [0.1, 0.1, 0.1]),
                              (None,       [0.2])])
    def test_returns_right_answers_for_valid_table(self, col2_arr, expected):
        srf = opt_surf.SpectralSurface()
        srf.table.add_column(Column(name="col1", data=[0.8]*len(expected)))
        if col2_arr:
            srf.table.add_column(Column(name="col2", data=col2_arr))

        col3 = srf._compliment_array("col1", "col2")
        if sys.version_info.major >= 3:
            assert col3.data == pytest.approx(expected)
            assert len(col3.data) == len(expected)
        else:
            warnings.warn("Data equality isn't tested for 2.7")


class TestQuantify:
    @pytest.mark.parametrize("item, unit, expected",
                             [(1, "m", 1*u.m),
                              (1, u.m, 1*u.m),
                              ([1], u.m, [1] * u.m),
                              (1 * u.m, u.m, 1 * u.m),
                              ([1] * u.m, u.m, [1] * u.m)])
    def test_return_quantity_for_valid_data(self, item, unit, expected):
        assert opt_surf.quantify(item, unit) == expected


@pytest.mark.usefixtures("unity_flux")
class TestMakeEmissionFromArray:
    @pytest.mark.parametrize("emission_unit",
                             ["ph s-1 m-2 um-1",
                              "ph s-1 m-2 Hz-1",
                              "erg s-1 m-2 um-1",
                              "jansky"])
    def test_source_spectrum_returned_for_synphot_units(self, emission_unit,
                                                        unity_flux):
        meta = {"emission_unit": emission_unit}
        out = opt_surf.make_emission_from_array(*unity_flux, meta=meta)
        assert isinstance(out, SourceSpectrum)

    @pytest.mark.parametrize("emission_unit",
                             ["ph s-1 m-2 um-1 arcsec-2",
                              "jansky sr-1"])
    def test_source_spectrum_returned_for_solid_angles(self, emission_unit,
                                                       unity_flux):
        meta = {"emission_unit": emission_unit}
        out = opt_surf.make_emission_from_array(*unity_flux, meta=meta)
        assert isinstance(out, SourceSpectrum)

    @pytest.mark.parametrize("emission_unit",
                             ["ph s-1 m-2",
                              "ph s-1 m-2 bin-1",
                              "erg s-1 m-2 arcsec-2",
                              "erg s-1 m-2 bin-1 sr-2"])
    def test_source_spectrum_returned_for_bin_units(self, emission_unit,
                                                    unity_flux):
        meta = {"emission_unit": emission_unit}
        out = opt_surf.make_emission_from_array(*unity_flux, meta=meta)
        assert isinstance(out, SourceSpectrum)


class TestNormaliseBinnedFlux:
    def test_returns_correct_normalisation(self):



        pass

