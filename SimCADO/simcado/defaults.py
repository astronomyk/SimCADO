'''Create default noise and PSF cubes'''
import sys
import os
import inspect
__pkg_dir__ = os.path.dirname(inspect.getfile(inspect.currentframe()))

import multiprocessing as mp

import numpy as np
import astropy.io.fits as fits

try:
    import simcado.detector as fpa
    import simcado.input as ig
except ImportError:
    import detector as fpa
    import input as ig

__all__ = ["make_noise_cube", "make_poppy_cube"]


def _make_noise(i):
    '''Single-plane noise realisation with HxRG noise'''
    pca_file = os.path.join(__pkg_dir__, "data", "FPA_nirspec_pca0.fits")

    ng_h4rg = fpa.HXRGNoise(naxis1=4096, naxis2=4096, naxis3=i,
                            n_out=32, nroh=8, pca0_file=pca_file)
    return ng_h4rg.mknoise(o_file=None, rd_noise=5, pedestal=4,
                           c_pink=3, u_pink=1, acn=0.5)


def make_noise_cube(num_layers=25, filename="noise.fits"):
    '''Create a cube with HxRG noise realisations'''
    if __name__ == "__main__":
        pool = mp.Pool(processes=mp.cpu_count()-1)
        frames = pool.map(_make_noise, np.ones(num_layers))
        hdu = fits.HDUList([fits.PrimaryHDU(frames[0])] + \
                           [fits.ImageHDU(frames[i])
                            for i in range(1, num_layers)])
        hdu.writeto(filename, clobber=True, checksum=True)
    else:
        frames = _make_noise(num_layers)
        hdu = fits.HDUList([fits.PrimaryHDU(frames[0])] + \
                           [fits.ImageHDU(frames[i])
                            for i in range(1, num_layers)])
        hdu.writeto(filename, clobber=True, checksum=True)


def make_poppy_cube():
    '''Create a cube with wavelength-dependent E-ELT PSFs'''
    # it makes no differnce if I use 3 or 7 cores. Average time per PSF is 24
    # seconds regardless of core number if cpus > 2
    if __name__ == "__main__":
        ig.poppy_eelt_psf_cube(lam_bin_centers=np.arange(0.7, 2.51, 0.01),
                               filename="poppy_cube.fits",
                               size=1023, pix_res=0.001,
                               cpus=mp.cpu_count() // 2)

if __name__ == "__main__":

    if "psf" in sys.argv[1:]:
        make_poppy_cube()

    if "noise" in sys.argv[1:]:
        make_noise_cube()
