"""
Extra utility functions that aren't needed for everyday SimCADO
"""    
    
    import numpy as np
    from astropy.io import fits

def centre_psfs(fname, output=None, x_cen=511.5, y_cen=511.5, tile_width=32):
    """
    Centre the PSFs in a PSF Cube FITS file made from external PSFs and write to disk
    
    Parameters
    ----------
    fname : str
        input file name
    output : str
        output file name - if None, fname is clobbered
    x_cen, y_cen : float
        [pixel] the position of where the new centre of the PSF will be. Default (511.5, 511.5)
    tile_width : int
        [pixel] the width of the tile used to centre the PSF

    Examples
    --------
    :
        >>> from simcado import __data_dir__
        >>>
        >>> fname = __data_dir__+"/PSF_SCAO.fits"
        >>> output = "PSF_SCAO_centred.fits"
        >>> centre_psfs(fname, output)

    """
    
    from scipy.ndimage import shift
    from astropy.stats import sigma_clipped_stats
    from photutils import DAOStarFinder
    
    x0 = x_cen - tile_width // 2
    y0 = y_cen - tile_width // 2 
    w = tile_width
        
    xc_pix = tile_width // 2 + x_cen % 1
    yc_pix = tile_width // 2 + y_cen % 1
    print(xc_pix, yc_pix)
    fits_file = fits.open(fname)

    for hdu in fits_file:
        im = hdu.data

        iters=5
        for i in range(iters):
            mean, median, std = sigma_clipped_stats(im[y0:y0+w, x0:x0+w], sigma=3.0, iters=5) 
            daofind = DAOStarFinder(fwhm=3.0, threshold=20.*std)    
            sources = daofind(im[y0:y0+w, x0:x0+w])

            i = np.argmax(sources["flux"])
            dx, dy = sources["xcentroid"][i] - xc_pix, sources["ycentroid"][i] - yc_pix

            im = shift(im, (-dy, -dx), order=1)
            im / np.sum(im)
            hdu.data = im
            
        print(sources)
    
    if output is None:
        fits_file.writeto(fname, clobber=True)
    else:
        fits_file.writeto(output, clobber=True)