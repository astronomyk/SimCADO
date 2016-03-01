"""InputGenerator"""
################################################
# InputGenerator
#
# Put various functions here that generate input data for a simulation
#
#
#



def analytic_eelt_psf_cube(lam_bin_centers, filename=None, **kwargs):
    """
    Generate a FITS file with E-ELT PSFs for a range of wavelengths.
    
    Parameters:
    ===========
    lam_bin_centers: [um] the centres of each wavelength bin
    filename: [None] path name to where the FITS file should be saved. If not 
              given, the file is stored in the working directory
    
    Optional Parameters:
    ====================
    pix_res: [arcsec] the angular resolution of the pixels. Default is 1 mas
    diameter: [m]
    """
    params = {  "diameter_out"  :37,
                "diameter_in"   :5.5,
                "flattoflat"    :1.45,
                "gap"           :0.004,
                "n_spiders"     :6,
                "pix_res"       :0.001,
                "size"          :255
                "oversample"    :1     }
    params.update(kwargs)
    
    rings = int(0.65 * params["diameter_out"] / params["flattoflat"])
    m1 = poppy.MultiHexagonAperture(rings=rings, 
                                    flattoflat=params["flattoflat"], 
                                    gap=params["gap"])
    pri = poppy.CircularAperture(radius=params["diameter_out"]/2)
    sec = poppy.SecondaryObscuration(secondary_radius=params["diameter_out"]/2,
                                     n_supports=params["n_spiders"], 
                                     support_width=0.5)
    eelt = poppy.CompoundAnalyticOptic( opticslist=[m1, pri, sec], name='E-ELT')

    osys = poppy.OpticalSystem()
    osys.addPupil(eelt)
    osys.addDetector(pixelscale=params["pix_res"],
                     fov_arcsec=params["pix_res"] * params["size"], 
                     oversample=params["oversample"])

    psfHDU = [osys.calcPSF(lam * 1E-6)[0] for lam in lam_bin_centers]
    hdulist = fits.HDUList(psfHDU)
    
    if filename is None:
        return hdulist
    else: 
        hdulist.writeto(filename, clobber=True)