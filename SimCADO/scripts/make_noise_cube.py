import Detector as fpa
import InputGenerator as ig
import multiprocessing as mp
import astropy.io.fits as fits

def make_noise(i):
    ng_h4rg = fpa.HXRGNoise(naxis1 = 4096, naxis2 = 4096, n_out = 32, nroh = 8, 
                            pca0_file = "../data/FPA_nirspec_pca0.fits")
    return ng_h4rg.mknoise(o_file = None, rd_noise = 5, pedestal = 4, 
                           c_pink = 3, u_pink = 1, acn = 0.5)


def make_noise_cube(n=None):                           
    
    if __name__ == "__main__":
        if n is None: 
            n=25
        pool = mp.Pool(processes=mp.cpu_count()-1)
        frames = pool.map(make_noise, range(n))
        hdu = fits.HDUList([fits.PrimaryHDU(frames[0])] + \
                           [fits.ImageHDU(frames[i]) for i in range(1,n)])
        hdu.writeto("noise.fits", clobber=True)
                           
                           
                           
def make_poppy_cube(lam, filename)
    if __name__ == "__main__":                       
                           


    
    
