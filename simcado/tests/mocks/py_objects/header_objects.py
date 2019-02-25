from astropy import wcs


def _basic_fov_header():
    w, h = 150, 150
    skywcs = wcs.WCS(naxis=2)
    skywcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    skywcs.wcs.cdelt = [0.1, 0.1]
    skywcs.wcs.cunit = ["arcsec", "arcsec"]
    skywcs.wcs.crval = [0, 0]
    skywcs.wcs.crpix = [w / 2, h / 2]

    detwcs = wcs.WCS(naxis=2, key="D")
    detwcs.wcs.ctype = ["LINEAR", "LINEAR"]
    detwcs.wcs.cdelt = [0.1, 0.1]
    detwcs.wcs.cunit = ["mm", "mm"]
    detwcs.wcs.crval = [0, 0]
    detwcs.wcs.crpix = [w / 2, h / 2]

    skyhdr = skywcs.to_header()
    dethdr = detwcs.to_header()
    skyhdr.update(dethdr)
    skyhdr["NAXIS"] = 2
    skyhdr["NAXIS1"] = w
    skyhdr["NAXIS2"] = h

    return skyhdr
