MAORY PSFs
==========

Here you find links to the October 2017 release of MAORY PSFs for both the SCAO
and MCAO modes, formatted for use by the new SimCADO v1.0 package.

.. warning:: MAORY PSFs are not in the public domain

    MAORY has only released these PSFs for internal use by the MICADO
    consortium. If you wish to publish anything using these PSFs, please
    **contact the MAORY team** regarding their publication policy.

.. note:: These PSFs are NOT formatted for legacy SimCADO (v0.5).

    These were put together for use with the **new SimCADO v1.0** code.
    Below you will find code to extracting the PSF kernels to be used with
    SimCADO v0.5.


MCAO PSFs
---------

The MCAO PSFs provide monochromatic PSFs kernels for IJHK filters at
(0.86, 1.245, 1.635, 2.145) micron and for 2 positions in the field: Near a
guide star (LGS or NGS), and far from a guide star.

* Wide-field (4mas / pixel) mode

  `<http://www.univie.ac.at/simcado/InstPkgSvr/psfs/MAORY_MCAO_FVPSF_4mas_20181203.fits>`_

* Zoom (1.5mas / pixel) mode

  `<http://www.univie.ac.at/simcado/InstPkgSvr/psfs/MAORY_MCAO_FVPSF_1.5mas_20181203.fits>`_


SCAO PSFs
---------

The SCAO PSFs provide monochromatic PSFs kernels for IzJHK filters at
(0.86, 1.02, 1.245, 1.635, 2.145) micron and for 49 positions in the field. The
positions are a 7x7 grid centred on the optical axis (0, 0), and spaced
7.5 arcsec in each direction - i.e. at (-22.5, -15, -7.5, 0, 7.5, 15, 22.5)
arcsec in each direction.

Each extension is a FITS image cube corresponds to one of the wavelengths noted
above, with each layer in the cube corresponding to PSF at the coordinates on
the field. The NGS is assumed to be on-axis.

PSF kernels are found in extensions 2 to 6. Extension 1 contains a table of the
coordinates with the corresponding cube layer index for the PSF kernel which
belongs at this coordinate.

* Wide-field (4mas / pixel) mode

  `<http://www.univie.ac.at/simcado/InstPkgSvr/psfs/MAORY_SCAO_FVPSF_4mas_20181203.fits>`_

* Zoom (1.5mas / pixel) mode

  `<http://www.univie.ac.at/simcado/InstPkgSvr/psfs/MAORY_SCAO_FVPSF_1.5mas_20181203.fits>`_


Using these PSFs with legacy SimCADO (v0.5)
-------------------------------------------
