Controlling simulations
=======================

This page describes two ways of simulating: the quick and the
repeatable. More about how SimCADO deals with commands can be found in
the `section on the UserCommands object <the-usercommands-object>`__.

Simulating the quick way
------------------------

The function ``sim.run()`` runs a simulation with all the default
parameters. That said, SimCADO has a safety switch to stop you wasting
10 minutes on a ``Source`` that will be invisible, or so bright that it
saturates the detector completely.

To get the full detector array, we set ``detector_array="full"``

::

    >>> im = sim.run(src, detector_array="full")

\*\* Note, by patient! \*\* Simulating 9x 4k detectors is heavy lifting.
Remeber to use the internals if you don’t plan on changing the optical
train.

Notes on simcado.run()
~~~~~~~~~~~~~~~~~~~~~~

There are two parameters still to be covered in detail. Please see the
documention (or read the docstring by using ``sim.run?``). They are:

-  ``filename=`` : instead of returning the FITS object to the console,
   save it to disk under this name
-  ``mode=`` : default is “wide” for 4mas imaging. Also accepts “zoom”
   for 1.5mas imaging mode.

Viewing a full detector read out
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each chip read-out is stored in a FITS extension. There are 3 options to
view the images:

1. Save the FITS object to disk and use DS9, etc

   ::

       >>> im.writeto("full_detector.fits")

2. Use the ``filename=`` in ``simcado.run()`` to save directly to disk,
   and bypass the console. This is useful for scripting with SimCADO.

   ::

       >>> simcado.run(src, filename="full_detector.fits")

3. Use the SimCADO function ``.plot_detector()`` from the ``.detector``
   module

   ::

       >>> im, (cmd, opt, fpa) = simcado.run(src, filename="full_detector.fits", 
                                               return_internals=True)
       >>> simcado.detector.plot_detector(fpa)   

The 3rd option is probably the least favourable as there are no options
available, but it allows you to see what the readout will look like in a
mosaic mode.

Using KEYWORD=VALUE pairs from the configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SimCADO is controlled with a series of keyword-value pairs contained in
a configuration file. **The defaults are the best approximation to
MICADO** so changing them is *not* recommended if you want to simulate
*MICADO* images. There are however some which are useful to play around
with.

!!! Note

::

    SimCADO is CAsE sEnsiTIve**. All KEYWORDS are writen with capital letters.

Similar to SExtractor, SimCADO provides a way to dump both commonly used
and less-common keywords to a file with command
``sim.commands.dump_defaults()``:

::

    >>> sim.commands.dump_defaults()
    #########################################################
    ##            FREQUENTLY USED KEYWORDS                 ##
    #########################################################

    #-------------------------------------------------------#
    # To dump a file with all keywords, use:                #
    # >>>  sim.commands.dump_defaults(filename, type="all") #
    #-------------------------------------------------------#

    OBS_EXPTIME             60          # [sec] simulated exposure time
    OBS_NDIT                1           # [#] number of exposures taken
    INST_FILTER_TC          Ks          # [<filename>, string(filter name)] List acceptable filters with >>> simcado.optics.get_filter_set()

    ATMO_USE_ATMO_BG        yes         # [yes/no]
    SCOPE_USE_MIRROR_BG     yes         # [yes/no]
    INST_USE_AO_MIRROR_BG   yes         # [yes/no]
    FPA_USE_NOISE           yes         # [yes/no]

    FPA_CHIP_LAYOUT         small       # [small/centre/full] description of the chip layout on the detector array. 
    SCOPE_PSF_FILE          scao        # ["scao" (default), <filename>, "ltao", "mcao", "poppy"] import a PSF from a file. Default is <pkg_dir>/data/PSF_SCAO.fits

    SIM_DETECTOR_PIX_SCALE  0.004       # [arcsec] plate scale of the detector
    SIM_VERBOSE             yes         # [yes/no] print information on the simulation run

    OBS_ZENITH_DIST         60          # [deg] from zenith
    INST_ADC_PERFORMANCE    100         # [%] how well the ADC does its job

To list all keyword-value pairs, use:

::

    >>>  sim.commands.dump_defaults(filename, type="all")

Any of the KEYWORD=VALUE pairs can be passed as extras to ``sim.run()``.
For example if we wanted to observe in J-band for 60 minutes, we would
pass:

::

    src = sim.source.source_1E4_Msun_cluster()
    im = sim.run(src, OBS_EXPTIME=3600, INST_FILTER_TC="J")

The jupyter notebook `my\_first\_sim.ipynb <my_first_sim.ipynb>`__ has
more exmples of this.

The UserCommands object
-----------------------

Behind the scenes of the ``simcado.run()`` command, three objects are
created:

-  a ``UserCommands`` object - for holding all the information on how a
   simulation should be run
-  an ``OpticalTrain`` object - which contains the models to describe
   each effect that needs to be simulated
-  a ``Detector`` object - commonly referred to as an ``fpa`` or Focal
   Plane Array. It describes the layout of the detectros and holds the
   observed images.

The ``UserCommands`` object is arguably the most important of these
three, because the other two need the keyword-value pairs contained
within the ``UserCommands`` object to correctly describe the optical
train and detector for the simulation.

A ``UserCommands`` object is created by reading in the defaults conifg
file (``defaults.config``) and then updating any of the keywords that
the user (or function) provides. For example, we can see all the default
keyword-value pairs by calling:

::

    >>> cmd = sim.UserCommands()

The ``UserCommands`` object contains 7 ordered dictionaries, one for
each topic and one general dictionary. Each can be referenced
individually, however all are updated when a value changes.

1. cmd.cmds - contains all keyword-value pairs
2. cmd.atmo - keyword-value pairs for the atmosphere
3. cmd.scope - keyword-value pairs for the telescope
4. cmd.inst - keyword-value pairs for the instrument (plus AO system)
5. cmd.fpa - keyword-value pairs for the dector array
6. cmd.obs - keyword-value pairs for the observation
7. cmd.sim - keyword-value pairs for the simulation

A ``UserCommands`` object can be used as a dictionary itself, although
technically all that happens is that it references the general
dictionary ``cmd.cmds``. For example

::

    >>> cmd["OBS_EXPTIME"] = 60

is exactly the same as either of the following two expressions

::

    >>> cmd.cmds["OBS_EXPTIME"] = 60
    >>> cmd.obs["OBS_EXPTIME"] = 60

Therefore for the sake of ease, we recommoned treating the
``UserCommands`` object as a dictionary and just using the default
syntax: ``cmd["..."] = xxx``

Saving and loading a ``UserCommands`` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Saving
^^^^^^

In case you have made changes to the values in a ``UserCommands`` object
that you would like to keep for next time, a ``UserCommands`` object can
be saved to disk with the following command:

::

    >>> cmd = sim.UserCommands()
    >>> cmd.writeto(filename="my_cmds.txt")

SimCADO writes out the dictionary in ASCII format.

Loading
^^^^^^^

Creating a ``UserCommands`` object based on a text file is as simple as
passing the file path:

::

    >>> cmd = sim.UserCommands("my_cmds.txt")

Special attributes
~~~~~~~~~~~~~~~~~~

A ``UserCommands`` object not only contains a dictionary of
keyword-value pairs, but also a select number of parameters pertaining
to the optical train for quick access. These include values for:

-  the exposure time for simulations: ``cmd.exptime``
-  the primary mirror: ``cmd.area``, ``cmd.diameter``
-  the wavelength vector for purely spectral data (i.e. transmission
   curves): ``cmd.lam``
-  the wavelength centres and edges for each spectral bin:
   ``cmd.lam_bin_centres``, ``cmd.lam_bin_edges``
-  the mirror configuration: ``cmd.mirrors_telescope``,
   ``cmd.mirrors_ao``
-  the detector plate scale and internal sampling resolutions:
   ``cmd.fpa_res``, ``cmd.pix_res``

Mirror and Detector configuration files
---------------------------------------

A quick note on the other files that SimCADO uses when creating an
optical train and the appropriate keywords

The detector array
~~~~~~~~~~~~~~~~~~

The detector array is described by a text file containing information on
the plate scale and the positions of the detector chips:

::

    >>> sim.commands.dump_chip_layout(path=None)
    #  id    x_cen    y_cen   x_len   y_len
    #        arcsec   arcsec   pixel  pixel
        4        0        0    4096    4096
        0  -17.084  -17.084    4096    4096
        1        0  -17.084    4096    4096
        2   17.084  -17.084    4096    4096
        3  -21.484        0    4096    4096
        5   17.084        0    4096    4096
        6  -17.084   17.084    4096    4096
        7        0   17.084    4096    4096
        8   17.084   17.084    4096    4096
        

This small file can be saved to disk by passing a filename to the
``path=`` parameters

::

    >>> sim.commands.dump_chip_layout(path="my_fpa.txt")

Any detector array can be provided to SimCADO, as long as the text file
follows this format. For example the HAWK-I detector array (4x
HAWAII-2RG) would look like this:

::

    #  id    x_cen    y_cen    x_len   y_len
    #        arcsec   arcsec   pixel   pixel
        0      -116     -116    2048    2048
        1       116     -116    2048    2048
        2      -116      116    2048    2048
        3       116      116    2048    2048
        

To pass a detector array description to SimCADO, use the
``FPA_CHIP_LAYOUT`` keyword:

::

    >>> cmd = sim.UserCommands()
    >>> cmd["FPA_CHIP_LAYOUT"] = "hawki_chip_layout.txt"

or pass is directly to the ``sim.run()`` command:

::

    >>> sim.run(... , FPA_CHIP_LAYOUT="hawki_chip_layout.txt", ...)

The mirror configurations
~~~~~~~~~~~~~~~~~~~~~~~~~

The mirror configuration can be dumped either to the screen or to disk
by using:

::

    >>> dump_mirror_config(path=None, what="scope")
    #Mirror     Outer   Inner   Temp
    M1          37.3    11.1    0.  
    M2          4.2     0.545   0.  
    M3          3.8     0.14    0.  
    M4          2.4     0.      0.  
    M5          2.4     0.      0.  

If ``path=None`` the contents of the default file are printed to the
screen. The parameter ``what`` is for the section of the optical train
that should be shown - either ``scope`` for the telescope, or ``ao`` of
the AO system. For most existing telescope, this parameter is
irrelevant. For the MICADO/MAORY setup however another six optical
surfaces are introduced into the system.

It is possible to specifiy different mirror configurations using a text
file with the same format as above. For example the VLT unit telescope
mirror config files would look like this:

::

    #Mirror Outer   Inner   Temp
    M1      8.2     1.0     0.
    M2      1.116   0.05    0.
    M3      1.0     0.      0.

To use this mirro config file in SimCADO use the keywords
``SCOPE_MIRROR_LIST`` and ``INST_MIRROR_AO_LIST``

::

    >>> cmd = sim.UserCommands()
    >>> cmd["SCOPE_MIRROR_LIST"] = "vlt_mirrors.txt"
    >>> cmd["INST_MIRROR_AO_LIST"] = "none"
