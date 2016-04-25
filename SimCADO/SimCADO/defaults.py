import sys
import numpy as np

import multiprocessing as mp
import astropy.io.fits as fits

import detector as fpa
import opticsGenerator as ig

def _make_noise(i):
    ng_h4rg = fpa.HXRGNoise(naxis1 = 4096, naxis2 = 4096, n_out = 32, nroh = 8,
                            pca0_file = "../data/FPA_nirspec_pca0.fits")
    return ng_h4rg.mknoise(o_file = None, rd_noise = 5, pedestal = 4,
                           c_pink = 3, u_pink = 1, acn = 0.5)


def make_noise_cube():
    if __name__ == "__main__":
        n=25
        pool = mp.Pool(processes=mp.cpu_count()-1)
        frames = pool.map(_make_noise, range(n))
        hdu = fits.HDUList([fits.PrimaryHDU(frames[0])] + \
                           [fits.ImageHDU(frames[i]) for i in range(1,n)])
        hdu.writeto("noise.fits", clobber=True)



def make_poppy_cube():
    # it makes no differnce if I use 3 or 7 cores. Average time per PSF is 24 
    # seconds regardless of core number if cpus > 2
    if __name__ == "__main__":
        ig.poppy_eelt_psf_cube(lam_bin_centers=np.arange(0.7,2.51,0.01), 
                               filename="poppy_cube.fits", 
                               size=1023, pix_res=0.001, cpus=mp.cpu_count() // 2)

if __name__ == "__main__":                               
                               
    if "psf" in sys.argv[1:]:
        make_poppy_cube()

    if "noise" in sys.argv[1:]:
        make_noise()