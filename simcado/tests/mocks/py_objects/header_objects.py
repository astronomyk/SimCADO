from astropy import wcs


def _basic_fov_header():
    w, h = 150, 150
    fovwcs = wcs.WCS(naxis=2)
    fovwcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    fovwcs.wcs.cdelt = [0.1, 0.1]
    fovwcs.wcs.cunit = ["arcsec", "arcsec"]
    fovwcs.wcs.crval = [0, 0]
    fovwcs.wcs.crpix = [w / 2, h / 2]

    fovhdr = fovwcs.to_header()
    fovhdr["NAXIS"] = 2
    fovhdr["NAXIS1"] = w
    fovhdr["NAXIS2"] = h

    return fovhdr
