Custom instrument packages for SimCADO
======================================

.. warning:: This format will only work with SimCADO v0.4, v0.5, and v0.6

    This way of specifying files and parameters has been depreciated for newer
    versions of SimCADO. A new API description will be available by June 2019.


SimCADO is a (reasonably) flexible framework for simulating telescope and
instrument optical trains. It was primarily developed for use with the ELT
telescope and the MICADO / MAORY instrument, however many other instruments
can be simulated with SimCADO if the proper configuration files are provided.

This page aims to give a short **overview of which files** SimCADO needs in
order to simulate a custom optical train.

.. note:: We're happy to help implement other instruments with SimCADO

    If you would like to use SimCADO to model another telescope / instrument
    system, and would like some help with how to do this, please contact
    the SimCADO team. We're more that happy to help out.


Config Files
------------
The default SimCADO configuration file can be found in the data folder of the
installation directory: ``simcado/data/default.config``. This folder can be
found by calling::

    import simcado
    simcado.__data_dir__

To create a custom instrument package, we will need to create our own
``.config`` file. We can dump the default file by calling::

    simcado.commands.dump_defaults("my_new_instrument.config", selection="all")

Next we will need to alter the values to fit our instrument. For some keywords
we will need to provide file paths to a (correctly formatted) file describing
the given effect. The most important files are described below.


File paths
----------

.. note:: Use ``simcado.__search_path__`` to store files in multiple places

    SimCADO contains a function (``find_file``) which looks for file names in
    the directories contained in the list ``simcado.__search_path__``.
    By default this list contains the working directory, and the two main
    folders in the installation directory::

        >>> simcado.__search_path__
        ['./', '<pkg_dir>/simcado', '<pkg_dir>/simcado/data']

    By adding directory names to the ``__search_path__`` list at the beginning
    of a python script/session, we can tell SimCADO to look in other directories
    for the relevant files.


Required Files
--------------

PSF files
+++++++++

Keywords : SCOPE_PSF_FILE

This file should contain a (series of) PSF kernel(s) as FITS image extensions.
Each extension must contain the keyword ``WAVE0`` which tells SimCADO the
wavelength for which the PSF kernel is relevant (in ``m`` or in ``um``).

For MICADO, the default PSF FITS file contains extensions for each of the major
broadband filters: I [0.87um], z [1.05um], J [1.25um], H [1.65um], Ks [2.15um].
For filters with central wavelengths which fall between these PSF kernels,
SimCADO takes the extension with the closest wavelength.

In case we need to include field-variations in the PSF cube, the new PSF format
is also accepted by SimCADO. It's description can be found here_ :

.. _here : https://telescopy.readthedocs.io/en/latest/design/Effect_Descriptions.html#field-varying-psfs


Mirror transmission curves
++++++++++++++++++++++++++

Keywords : SCOPE_M1_TC, INST_MIRROR_AO_TC, INST_MIRROR_TC, INST_FILTER_TC

TC (transmission curve) files are ASCII table files with two columns. The
transmission profile can be as simple or as detailed as we want. Here is an
example of a very simple filter transmission file::

    # TC_filter_Ks.dat
    Wavelength  Transmission
    0.1         0
    1.89        0
    1.9         1
    2.4         1
    2.41        0
    3.0         0

SimCADO uses ``astropy.table`` to read the file, so if ``astropy`` can read it,
SimCADO can too.


Mirror lists
++++++++++++

Keywords : SCOPE_MIRROR_LIST, INST_MIRROR_AO_LIST

These list files contain information on the geometry of the mirrors included
in either the telescope section, or AO section of the optical system.
Below we have the first 3 lines of the ELT's mirror list file::

    # EC_mirrors_scope.tbl
    Mirror      Outer   Inner   Angle  Temp    Coating
    M1          37.3    11.1    0.     0.      TC_mirror_mgf2agal.dat

where Outer, Inner [m] refer to the diameter of the mirror, Angle [deg] is the
angle at which the plane of the mirror differs from that of the light beam,
Temp [deg C] is the mirror temperature for determining the thermal background,
and Coating is the filename of the Transmission Curve (TC) file describing
the spectral response of the mirror coating material (see above: Mirror
transmission curves).

The internal mirror configuration of the actual instrument is a special case.
It is assumed that the cryostat temperature is sufficiently low that the
blackbody emission from the internal mirrors is negligible. Hence only the
mirror transmission is important and can therefore be described with the
following two parameters::

    INST_NUM_MIRRORS
    INST_MIRROR_TC

For instruments without AO module, we can force SimCADO to ignore the AO
components by setting the following keyword::

    INST_USE_AO_MIRROR_BG = False


Wavefront errors
++++++++++++++++

Keywords : INST_WFE

The current version of SimCADO (v0.6) includes only the reduction in strehl
ratio induced by wave front error. It does this by calculating the peak of a
normalised 2D Gaussian kernel (i.e. between 1 and 0) for each wavelength in the
final system transmission curve, and applying an effective transmission loss
based on the Gaussian peak value. This is very much a 'first order
approximation' to including wavefront errors in the optical train. It should
be noted that this has been improved upon in SimCADO v1.0.

The file passed to ``FPA_QE`` should follow this table format (taken from the
MICADO+ELT ``INST_wfe.tbl`` file)::

    # wfe_rms   no_surfaces     material    optics
    # [nm]      [#]             [str]       [str]
    20          11              gold        mirror
    10          4               glass       entrance_window
    10          2               glass       filter
    10          8               glass       ADC


Optional files
--------------

Atmospheric spectra
+++++++++++++++++++

Keywords : ATMO_TC, ATMO_EC

By default SimCADO uses precalculated output from the ESO skycalc_ tool for the
atmospheric emission and transmission curves. The default tables
(``EC_sky_25.tbl`` and ``TC_sky_25.tbl``) are for PWV = 2.5mm and include
columns for various airmasses over a wavelength range or [0.3, 3.0]um.

.. _skycalc : https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC`

Additional transmission curves
++++++++++++++++++++++++++++++

Keywords : INST_DICHROIC_TC, INST_ENTR_WINDOW_TC, INST_PUPIL_TC

In case the transmission properties of any dichroics, entrance windows or pupil
optics need to be included. The file format is the same as above for
`Mirror transmission curves`_

Detector noise images
+++++++++++++++++++++

Keywords : FPA_NOISE_PATH

By default SimCADO uses a FITS file containing an image of the detector noise
pattern for the HAWAII-4RG detectors. We can include pre-calculated noise maps
by passing a FITS file to ``FPA_NOISE_PATH``. No special header keywords are
needed.

Detector linearity curve
++++++++++++++++++++++++

Keywords : FPA_LINEARITY_CURVE

Detector linearity can be included by passing an ASCII table with two columns
which relate the real incoming flux to the measured photon flux as measured by
the read-out electronics. The table should look like this::

    # real_flux measured_flux
    0           0
    1           1
    1000        998
    3500        3200
    ...         ...
    200000      180000
    1000000     180000

.. note:: Linearity is applied to any imaging observation, regardless of length

    Yes, this shouldn't be, but we haven't got around to fixing that yet. Hence
    to model long exposure observations (i.e. >1 min), it's best just to set
    ``FPA_LINEARITY_CURVE = "none"``

