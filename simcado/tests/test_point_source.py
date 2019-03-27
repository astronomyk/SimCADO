"""
Tests for:

source.point_source()
source.redshfit()
"""

import numpy as np
import simcado
import pytest


@pytest.mark.parametrize("sed", ["A0V", "spiral", "elliptical"])
@pytest.mark.parametrize("filter_name", ["J", "Ks"])
@pytest.mark.parametrize("magnitude", np.linspace(10, 20, 5))
def test_point_source_inputs(sed, filter_name, magnitude):
    """
    Testing that inputs with strings (filenames, SEDs) return the same values
    that inputs with Emission- and TransmissionCurves

    """

    filter_tc = simcado.optics.get_filter_curve(filter_name)
    lam, flux = simcado.source.SED(spec_type=sed, filter_name=filter_name, magnitude=magnitude)
    ec = simcado.source.scale_spectrum(lam, flux, mag=magnitude, filter_name=filter_name, return_ec=True)

    src1 = simcado.source.point_source(spectrum=sed, mag=magnitude, filter_name=filter_name)
    src2 = simcado.source.point_source(spectrum=ec, mag=magnitude, filter_name=filter_tc)

    nphot1 = src1.photons_in_range(lam_min=1.0, lam_max=2.4)
    nphot2 = src2.photons_in_range(lam_min=1.0, lam_max=2.4)

    assert np.isclose(nphot1, nphot2, 1e-3)


@pytest.mark.parametrize("spec_type", ["A0V", "G2V", "K0III"])
@pytest.mark.parametrize("filter_name", ["J", "Ks"])
@pytest.mark.parametrize("mag", np.linspace(10, 20, 5))
def test_vs_source_star(spec_type,  mag, filter_name):
    """
    Testing that source.star and source.point_source return the same answers
    Only works with stellar spectra ofc
    """
    src1 = simcado.source.star(spec_type, mag, filter_name)
    src2 = simcado.source.point_source(spec_type, mag, filter_name)
    nphot1 = src1.photons_in_range(lam_min=1.0, lam_max=2.4)
    nphot2 = src2.photons_in_range(lam_min=1.0, lam_max=2.4)

    assert np.isclose(nphot1, nphot2, 1e-3)


@pytest.mark.parametrize("z", [1, 2])
@pytest.mark.parametrize("spectrum", ["A0V", "elliptical"])
@pytest.mark.parametrize("mag", np.linspace(10, 20, 4))
@pytest.mark.parametrize("filter_name", ["J", "Ks"])
def test_redshift_inputs(z, spectrum, mag, filter_name):
    """
    Testing that inputs with strings (filenames, SEDs) return the same values
    that inputs with Emission- and TransmissionCurves

    """

    filter_tc = simcado.optics.get_filter_curve(filter_name)
    lam, flux = simcado.source.SED(spec_type=spectrum, filter_name=filter_name, magnitude=mag)
    ec = simcado.source.scale_spectrum(lam, flux, mag=mag, filter_name=filter_name, return_ec=True)

    ec1 = simcado.source.redshift_SED(z, spectrum, mag, filter_name)
    ec2 = simcado.source.redshift_SED(z, ec, mag, filter_tc)

    nphot1 = np.sum(ec1.val[(ec1.lam > 1) & (ec1.lam < 2.4)])
    nphot2 = np.sum(ec2.val[(ec1.lam > 1) & (ec2.lam < 2.4)])

    assert np.isclose(nphot1, nphot2, 1e-2)


