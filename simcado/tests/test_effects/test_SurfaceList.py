import os
import pytest
from pytest import approx

import numpy as np
from astropy import units as u

from synphot import SpectralElement, SourceSpectrum

import simcado as sim
from simcado.optics.effects import TERCurve
from simcado.optics.effects.surface_list import SurfaceList
from simcado.optics.radiometry import RadiometryTable
from simcado.optics.surface import SpectralSurface

from matplotlib import pyplot as plt

PLOTS = False

MOCK_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                         "../mocks/MICADO_SCAO_WIDE/"))
sim.rc.__search_path__ += [MOCK_PATH]


class TestSurfaceListInit:
    def test_initialise_with_nothing(self):
        assert isinstance(SurfaceList(), SurfaceList)

    def test_initalises_with_list_of_surfaces(self):
        filename = os.path.join(MOCK_PATH, "EC_mirrors_MICADO_Wide.tbl")
        surf_list = SurfaceList(filename=filename)
        assert isinstance(surf_list, SurfaceList)


class TestRadiometryTableAttribute:
    def test_returns_radiometry_table_object(self):
        filename = os.path.join(MOCK_PATH, "EC_mirrors_MICADO_Wide.tbl")
        etendue = 5776 * u.m ** 2 * u.mas ** 2
        surf_list = SurfaceList(filename=filename, etendue=etendue)

        assert isinstance(surf_list.radiometry_table, RadiometryTable)
        assert isinstance(surf_list.get_emission(), SourceSpectrum)

        if PLOTS:
            wave = np.arange(10, 200.51, 1) * u.um
            plt.plot(wave, surf_list.get_emission()(wave))
            plt.semilogy()
            plt.show()


class TestAddSurfaceList:
    def test_second_list_is_joined_to_first(self):
        etendue = 5776 * u.m ** 2 * u.mas ** 2

        filename1 = os.path.join(MOCK_PATH, "EC_mirrors_MICADO_Wide.tbl")
        surf_list1 = SurfaceList(filename=filename1, etendue=etendue)
        len1 = len(surf_list1.radiometry_table.table)

        filename2 = os.path.join(MOCK_PATH, "EC_mirrors_ELT.tbl")
        surf_list2 = SurfaceList(filename=filename2, etendue=etendue)
        len2 = len(surf_list2.radiometry_table.table)

        surf_list1.add_surface_list(surf_list2, prepend=True)
        len3 = len(surf_list1.radiometry_table.table)

        assert len3 == len2 + len1


class TestAddSurface:
    def test_extra_surface_is_joined_to_list(self):
        etendue = 5776 * u.m ** 2 * u.mas ** 2
        surf_list = SurfaceList(filename="EC_mirrors_MICADO_Wide.tbl",
                                etendue=etendue, name="MICADO Mirror List")
        surf = TERCurve(filename="TC_filter_Ks.dat", name="filter",
                        action="transmission", outer=0.1, temp=0)
        len1 = len(surf_list.radiometry_table.table)
        surf_list.add_surface(surf, "filter")

        assert len(surf_list.radiometry_table.table) == len1 + 1

        if PLOTS:
            wave = np.linspace(0.5, 2.5, 100)*u.um
            plt.subplot(121)
            plt.plot(wave, surf_list.get_throughput()(wave))

            plt.subplot(122)
            plt.plot(wave, surf_list.get_emission()(wave))
            plt.semilogy()
            plt.show()


class TestTERCurveInit:
    def test_initialise_with_nothing(self):
        assert isinstance(TERCurve(), TERCurve)

    def test_initalises_with_list_of_surfaces(self):
        filename = os.path.join(MOCK_PATH, "TC_filter_Ks.dat")
        surf = TERCurve(filename=filename)
        assert isinstance(surf, TERCurve)

    def test_initalises_with_two_arrays(self):
        surf = TERCurve(wavelength=np.array([0.5, 2.5]),
                        transmission=np.array([1, 1]),
                        wavelength_unit="um")
        assert isinstance(surf, TERCurve)


class TestSurfaceAttribute:
    def test_returns_surface_object(self):
        filename = os.path.join(MOCK_PATH, "TC_filter_Ks.dat")
        surf = TERCurve(filename=filename)

        assert isinstance(surf.surface, SpectralSurface)
        assert isinstance(surf.surface.transmission, SpectralElement)
        assert isinstance(surf.surface.emission, SourceSpectrum)

    def test_returns_surface_object_for_arrays(self):
        surf = TERCurve(wavelength=[0.5, 1.5, 2.5],
                        transmission=[0.1, 0.1, 0.1],
                        wavelength_unit="um")

        assert isinstance(surf.surface, SpectralSurface)
        assert isinstance(surf.surface.transmission, SpectralElement)
        assert isinstance(surf.surface.emission, SourceSpectrum)
