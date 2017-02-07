import numpy as np
import simcado as sim

def test_stars_params():
    # Test combinations of input parameters
    assert sim.source.stars()
    assert sim.source.stars("G2V", 15)
    assert sim.source.stars("G2V", [15, 12])
    assert sim.source.stars(['G2V', 'F2V'], 15)
    assert sim.source.stars(["G2V", "F2V"], [15, 12])


def test_source_resample_equivalency():

    n=8

    im = np.ones((100,100))
    lam, spec = sim.source.SED("M0V", "K", 20)
    src = sim.source.source_from_image(im, lam, spec, 0.004, oversample=n)
    hdu, (cmd, opt, fpa) = sim.run(src, return_internals=True)

    im = np.ones((n*100, n*100)) / (n**2)
    lam, spec = sim.source.SED("M0V", "K", 20)
    src2 = sim.source.source_from_image(im, lam, spec, 0.004/n, oversample=1)
    hdu2, (cmd2, opt2, fpa2) = sim.run(src2, return_internals=True)

    diff = np.sum(np.abs(fpa2.chips[0].array) - np.abs(fpa.chips[0].array)) / \
                                                    np.sum(fpa.chips[0].array)

    assert diff < 1E-4


def test_stars_delivers_the_same_as_SED():

    #sim.source.stars([])

    spec_types = ["A0V", "A0V", "M5V"]
    lam, spec = sim.source.SED(spec_type=spec_types[0],
                               filter_name="Ks",
                               magnitude=20.)
    vega_SED  = sim.spectral.EmissionCurve(lam=lam, val=spec)
    vega_star = sim.source.stars(spec_types=spec_types[0],
                                 filter_name="Ks",
                                 mags=20)

    assert np.sum(vega_SED.val) == np.sum(vega_star.spectra[0] * \
                                                            vega_star.weight[0])
