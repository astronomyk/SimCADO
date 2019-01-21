import numpy as np
import astropy
import simcado
import astropy.units as u
import os


#### Helper functions to check against tabulated values
def create_vega_table():
    """
    Returns astropy.table with tabulated values for Vega star
    Data taken from http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
    -------

    """
    filters = astropy.table.Column(name="filters",
                                   data=["U", "B", "V",	"R", "I", "J", "H", "K"])
    lambda_eff = astropy.table.Column(name="lambda_eff",
                                      data=[0.36, 0.438, 0.545, 0.641, 0.798, 1.22, 1.63, 2.19],
                                      unit=u.um)
    delta_lambda = astropy.table.Column(name="delta_lambda",
                                        data=[0.06, 0.09, 0.085, 0.15, 0.15, 0.26, 0.29, 0.41],
                                        unit=u.um)
    f_nu = astropy.table.Column(name="f_nu",
                                data=1e-20*np.array([0.79, 4.063, 3.636, 3.064, 2.416, 1.589, 1.021, 0.64]),
                                unit=u.Unit("erg / (s cm2 Hz)"))
    f_lambda = astropy.table.Column(name="f_lambda",
                                    data = 1e-11*np.array([417.5, 632, 363.1, 217.7, 112.6, 31.47, 11.38, 3.961]),
                                    unit=u.Unit("erg / (s cm2 angstrom)"))
    n_ph = astropy.table.Column(name="n_photons",
                                data=[756.1, 1392.6, 995.5, 702.0, 452.0, 193.1, 93.3, 43.6],
                                unit=u.Unit("photon / (s cm2 angstrom)"))
    vega_table = astropy.table.Table((filters,lambda_eff,delta_lambda,f_nu,f_lambda,n_ph),
                                     meta={"source" :"http://www.astronomy.ohio-state.edu/~martini/usefuldata.html"})
    return vega_table


def test_returning_numbers(filter_name="TC_filter_K.dat"):
    """
    Check that zero_magnitude_photon returns a positive number
    Parameters
    ----------
    filter_name

    Returns
    -------

    """

    filter_file = os.path.join("mocks", filter_name)
    nph = simcado.source.zero_magnitude_photon_flux(filter_file)
    assert nph > 0


def test_nph_from_sources(filter_name="TC_filter_K.dat"):
    """
    Check that differences are below 1% level if a file or a transmission curve is passed

    Parameters
    ----------
    filter_name

    Returns
    -------

    """

    filter_file = os.path.join("mocks", filter_name)
    nph_from_filter = simcado.source.zero_magnitude_photon_flux(filter_file)
    tc = simcado.optics.get_filter_curve(filter_name)
    nph_from_tc = simcado.source.zero_magnitude_photon_flux(tc)
    assert np.abs((nph_from_filter / nph_from_tc) - 1) < 0.01


def test_mag_to_photons_to_mag(filter_name="TC_filter_K.dat", magnitude=0):
    """
    Test that photons_to_mag and mag_to_photons return the original value
    (as they use zero_magnitude_photon_flux for their calculations)
    Parameters
    ----------
    filter_name
    magnitude

    Returns
    -------
    """

    filter_file = os.path.join("mocks", filter_name)
    nph = simcado.source.mag_to_photons(filter_file, magnitude)
    mag = simcado.source.photons_to_mag(filter_file, nph)
    assert magnitude == mag

def test_create_vega_table():
    table = create_vega_table()
    assert isinstance(table,astropy.table.Table)

def test_photon_flux_in_filter(filter_name="K"):
    """
    Compare tabulated photon fluxes with those returned by simcado.source.zero_magnitude_photon_flux
    Tolerance 10% to account for different filter curves shapes

    Parameters
    ----------
    filter_name: choices are U, B, V, R, I, J, H, K

    Returns
    -------

    """
    filter_file = os.path.join("mocks", "TC_filter_" + filter_name + ".dat")
    n_photons_simcado = simcado.source.zero_magnitude_photon_flux(filter_file)
    tabulated_fluxes = create_vega_table()
    mask = tabulated_fluxes["filters"] == filter_name
    t_temp = tabulated_fluxes[mask]
    n_photons_tab = t_temp["n_photons"]
    delta_lambda_tab = t_temp["delta_lambda"]
    n_photons_tab = n_photons_tab.to(u.Unit("photon / (um m2 s)")) * delta_lambda_tab
    assert np.abs((n_photons_simcado / n_photons_tab[0].value) - 1) < 0.1


# test_returning_numbers("TC_filter_Pa-beta.dat")
# test_nph_from_sources("TC_filter_H.dat")
# test_mag_to_photons_to_mag("TC_filter_H.dat", -10)
# test_create_vega_table()
# test_photon_flux_in_filter("I")
