Updates to SimCADO
==================

Data updates
------------

The data files used by SimCADO are continually being updated to reflect
the newest information coming from the MICADO work packages.

To update your version of SimCADO, use the ``simcado.get_extras()``
command.

2019-02-19
~~~~~~~~~~

- Documentation Updated
- Reference spectrum now based on synphot


2018-11-23
~~~~~~~~~~

- Moved the documentation to ReadTheDocs (finally)
- Added filter plotting functionality


2017-01-31
~~~~~~~~~~

**New stable version available: SimCADO 0.4**

2017-01-18
~~~~~~~~~~

-  Added E-ELT mirror transmission curve TC\_mirror\_EELT.dat

2017-01-05
~~~~~~~~~~

-  VLT mirror coating TC\_aluminium.dat
-  HAWK-I filter curves: TC\_filter\_J.dat
   TC\_filter\_H.dat
   TC\_filter\_Ks.dat
   TC\_filter\_Y.dat
-  Updated detector layout FPA\_chip\_layout.dat

2016-11-12
~~~~~~~~~~

-  ``TC_surface.dat`` Updated with more optimistic Strehl values for JHK
   central wavelengths
-  ``TC_mirror_gold.dat`` Added a reflectivity curve for an unprotected
   gold coating
-  ``default.config`` Changed ``INST_MIRROR_TC`` to reference
   ``TC_mirror_gold.dat``

Source Code updates
-------------------

The current stable version is SimCADO v0.4. The current development
version is SimCADO v0.5dev

Development version
~~~~~~~~~~~~~~~~~~~

2017-01-05 \* SimCADO now preloads the transmission curves. Generating
an OpticalTrain now only takes ~2 seconds (compared to 20 previously)

2016-12-06 \* Added functionality to generate “ideal” AO PSFs using
POPPY for the diffraction limited core and added a Seeing halo

2016-11-12 \* Bug fix in simcado.psf - psf.lam\_bin\_centers was taking
opt\_train.lam\_bin\_centers instead of using the “WAVE0” keywords in
the FITS header \* Bug fix in source.psf - apply\_optical\_train() was
asking for PSFs outisde of the opt.train.psf range. related to the point
above
