Getting Started with SimCADO
============================

SimCADO can be super easy to use, or super complicated. The level of
complexity is completely up to the user. A basic simulation involves
only 1 thing: a :class:`.Source` object to describe the observable object.
Once the user has created this object, the function 
:func:`.simcado.simulation.run()` is all that needs to be called. Controlling 
the parameters of the simulation can be done either by passing keyword-value 
pairs, or my using a :class:`.UserCommands` dictionary.

.. contents::


.. note:: 
    Even though this documentation is not yet complete (it is a very big job), 
    lots **more information is included in the docstrings** of every SimCADO 
    function and class. These can be easily viewed in the interactive python 
    interface (iPython, or Jupyter Notebook) with either the question mark 
    operator or by using ``SHIFT+TAB`` with the cursor over the function name.

Doing this in iPython will call up the docstring: ::

    >>> simcado.Source?

Source
------

For full details, please see the :doc:`API <simcado>` and examples of
:doc:`Source Objects <C_Source>`

The :class:`~.Source` class is probably the most important class for testing
science cases. Therefore spending time on creating accurate :class:`~.Source`
representations of the object of interest is key to getting good results
with SimCADO. :class:`~.Source` objects can be created from scratch, with
functions provided by SimCADO, or by loading in a pre-existing
:class:`~.Source`-FITS file.

.. note:: 
    **SimCADO is CaSe SensITIVe!** SimCADO has the class ``simcado.Source()``
    and the module :mod:`simcado.source`. These should not be confused. 
    :mod:`simcado.source` is the module which contains the class :class:`.Source` 
    and all the helper functions for creating various types of :class:`.Source` 
    objects. The source code for the  class :class:`.Source` is actually in
    :class:`simcado.source.Source`, however to make things easy, :class:`.Source`
    is available directly as ``simcado.Source()``. Be careful and remember 
    ``simcado.Source != simcado.source``.

For a description of the :class:`~.Source` object, and the :mod:`.source` module,
see :doc:`How SimCADO works <A_DeepStuff#Source>`.

Loading a pre-existing :class:`.Source` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To load in a pre-existing :class:`.Source` (i.e. one that you saved earlier),
specify the keyword ``filename=`` when initialising the :class:`.Source`
object.::

    >>> import simcado as sim
    >>> my_src = sim.Source(filename="star_grid.fits")

:class:`.Source`-FITS files have a very specific file format, so it’s best to
only import files that were generated directly from other :class:`.Source`
objects. It’s a chicken/egg scenario, which is why the next section
deals with creating :class:`.Source` objects in memory. For a description of
the file format for saved :class:`.Source` objects, see :ref:`“File Format of saved
Source objects” <deepstuff>`.


Making a :class:`.Source` with SimCADO’s in-built functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :mod:`simcado.source` module provides an ever-increasing series of
functions to create :class:`.Source` objects in memory. These include, (from
:mod:`simcado.source`)

-  :func:`empty_sky() <.empty_sky>`
-  :func:`star(mag, filter_name="K", ...) <.star>`
-  :func:`stars(mags, x, y, ...) <.stars>`
-  :func:`cluster(mass=1E4, distance=50000, ...) <.cluster>`
-  :func:`source_from_image(images, lam, spectra, pix_res, ...) <.source_from_image>`

Two useful functions here are :func:`.stars()` and :func:`.source_from_image()`

-  :func:`.stars` takes a list of magnitudes (and optionally spectral
   types) and positions for a common broad-band filter (default is “K”)
   and generates a :class:`.Source` object with those stars in the field. ::

        >>> x, y = [-2.5, 0.7, 16.3], [3.3, -0.2, 25.1]
        >>> mags, spec_types = [25,21,28], ["K0V", "A0III", "G2V"]
        >>> filt = "H"
        >>>
        >>> my_src = sim.source.stars(mags=mags, x=x, y=y, filter_name=filt, 
                                                        spec_types=spec_types)

-  :func:`.source_from_image` creates a :class:`.Source` based on a 2D numpy
   array provided by the user. The 2D array can come from anywhere,
   e.g. the data from a FITS image, a BITMAP image, from memory, etc.
   Alongside the image, the user must provide a spectrum (plus a vector
   with the bin centres) and the pixel field of view (e.g. 0.004 arcsec
   for MICADO). SimCADO then extracts all pixels from the image which
   have values above ``flux_threshold`` (defualt is 0) and saves these
   pixel coordinates. The spectrum provided is then connected to these
   pixel, and scaled by the pixel value. ::

        >>> # ... Create an image - a circle with a radius of 20 pixels on a 
        >>> # ... grid 200 pixel wide
        >>> XX = np.array([np.arange(-100,101)]*201) 
        >>> im = np.sqrt(XX**2 + XX.transpose()**2)
        >>> im[im>20] = 0; im[im>0] = 1
        >>>
        >>> # ... Pull in the spectrum for a G2V star with K=20
        >>> lam, spec = simcado.source.SED("G2V", filter_name="K", magnitude=20)
        >>>
        >>> # ... Make the source object
        >>> my_src = sim.source.source_from_image(images=im, lam=lam, spectra=spec, plate_scale=0.004)

SimCADO also provides a series of spectra for stars and galaxies,
however these are meant as a guide to those who are just starting out.
For serious work, the user is encouraged to provide their own spectra.
More information on the in-built spectra can be found in the
:doc:`Source Objects example <C_Source>` section.


Simulating with SimCADO
-----------------------

The quick, the dirty and the ugly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As seen on the :doc:`index <../index>` page, a simulation can be run using 3
lines of code: ::

    >>> import simcado
    >>> src = simcado.Source(filename="my_source.fits")
    >>> simcado.run(src, filename="my_image.fits")   

The :func:`.run` function is quite powerful. Many users may find that they
don’t need anything else to run the simulations they need. The full
function call looks like this: ::

    simcado.run(src, filename=None, 
                mode="wide", detector_layout="small",  
                cmds=None, opt_train=None, fpa=None, 
                return_internals=False,
                **kwargs)
               

Lets pull this function call apart in order of importance to the
simulation:

1. ``src``: Obviously the more important aspect is the :class:`.Source`
   object. Without a :class:`.Source` these is nothing to observe
2. ``filename``: Where to save the output FITS file. If ``None`` is
   provided (or the parameter is ignored), the output is returned to the
   user. This comes in handy if you are working in a
   Jupyter Notebook and wand to play with the output data
   immediately. Or if you are scripting with SimCADO and don’t want to
   be slowed down by writing all images to disk
3. Two important parameters here are ``mode`` and ``detector_layout``:
   These two define the MICADO observing modes.

Currently ``mode`` can be either ``"wide"`` (4mas/pixel) or ``"zoom"``
(1.5mas/pixel).

The ``detector_layout`` can also be changed to speed up simulations of
single objects. For example if the galaxy you’re interested in is at
z=5, you don’t need to read out all 9 MICADO chips for each observation.
In fact, a 1024x1024 window at the centre of the middle chip will
probably be enough. Therefore SimCADO offers the following “layouts” for
the detector - “small”, “wide”, “full”. The default is “small”.

-  ``small`` - 1x 1k-detector centred in the FoV
-  ``centre`` - 1x 4k-detector centred in the FoV
-  ``full`` - 9x 4k-detector as described by the keyword ``FPA_CHIP_LAYOUT``

1. ``cmds, opt_train, fpa`` are all parameters that allow you to provide
   custom built parts of the machinary. Say you have a set of commands
   saved from a previous simulation run which differ from the default
   values, then you can use these by passing a :class:`~.UserCommands` object
   via the ``cmd`` parameter. The same goes for passing an custom
   ``OpticalTrain`` object to ``opt_train`` and a custom ``Detector``
   object to ``fpa``. For more information see the relevant examples
   sections -
   :class:`~.UserCommands` :doc:`examples <B_UserCommands>`,
   :class:`~.OpticalTrain` :doc:`examples <D_OpticalTrain>`,
   :class:`~.Detector` :doc:`examples <E_Detector>`.

2. ``return_internals`` allows you to do the opposite of the previous
   three parameters. If you would like to save the :class:`~.UserCommands`,
   ``Detector`` and/or ``OpticalTrain`` from your simulation run, the by
   setting ``return_internals=True``, SimCADO will return these along
   with the simulated imagery. **Note** that this only works if
   ``filename=None``.

3. ``**kwargs``: Although ``kwargs`` is the last parameter, it actually
   allows you to control every aspect of the simulation. ``kwargs``
   takes any keyword-value pair that exist in the SimCADO configuration
   file, and so you can control single aspects of the simulation by
   passing these keyword-value pairs to :func:`.run`. For example, you can
   increase the exposure time of the image by passing ::

       simcado.run(src, … , OBS_EXPTIME=600, INST_FILTER_TC=“J”, …)

A list of all the available keyword-value pairs can be found in the
:doc:`Keywords section <Z_Keywords>` .

Alternatively you can dump a copy of the default parameters by calling
``simcado.commands.dump_defaults()``.


Changing Filters
~~~~~~~~~~~~~~~~

The keyword ``INST_FILTER_TC`` allows you to supply either the name of a
filter (i.e. “Ks”, “PaBeta”) or a path to an ASCII file containing a
filter curve. ``INST_FILTER_TC`` can be passed to :func:`.run` just like
any other SimCADO configuration keyword: ::

    >>> simcado.run(src, INST_FILTER_TC="J")
    >>> simcado.run(src, INST_FILTER_TC="path/to/my_filter_curve.txt")

SimCADO has some generic filters built in. These include all the regular
NIR broadband filters (I, z, Y, J, H, K, Ks). There are also some narrow
band filter. As the MICADO filter set is expected to change, we will not
list the SimCADO filter set here. Instead the user can find out which
filters are available by calling the function (as of Nov 2016): ::

    >>> print(sim.optics.get_filter_set())
    ['B', 'BrGamma', 'CH4_169', 'CH4_227', 'FeII_166', 'H', 'H2O_204', 'H2_212', 
     'Hcont_158', 'I', 'J', 'K', 'Ks', 'NH3_153', 'PaBeta', 'R', 'U', 'V', 'Y', 
     'z']

If you’d like to use your own filter curve, note that the ASCII file
should contain two columns - the first holds the wavelength values and
the second hold the transmission values between 0 and 1.

Setting the observation sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The important keywords here are: ``OBS_EXPTIME``, ``OBS_NDIT``

-  ``OBS_EXPTIME`` [in seconds] sets the length of a single exposure.
   The default setting is for a 60s exposure
-  ``OBS_NDIT`` sets how many exposures are taken. The default is 1.

Depending on what your intended use for SimCADO is, the keyword
``OBS_SAVE_ALL_FRAMES=["no", "yes"]`` could also be useful. The default
is to **not** save all the individual exposzures, but stack them and
return a single HDU object (or save to a single FITS file). If
``OBS_SAVE_ALL_FRAMES="yes"``, then a ``filename`` must also be given so
that each and every DIT can be saved to disk.

Reading out the detector
~~~~~~~~~~~~~~~~~~~~~~~~

**Warning**: running a full simulation could take ~10 minutes, depending
on how much RAM you have available::

    >>> simcado.run(detector_layout="small"")

The ``detector_layout`` keyword is key: ::

    detector_layout : str, optional
        ["small", "centre", "full"] Default is "small".

Where each of the strings means:

-  ``"small"`` - 1x 1k-detector centred in the FoV
-  ``"centre"`` - 1x 4k-detector centred in the FoV
-  ``"full"`` - 9x 4k-detector as per MICADO imaging mode (either 4mas
   or 1.5mas)
-  ``"default"`` - depends on “mode” keyword. Full MICADO 9 chip
   detector array for either 4mas or 1.5mas modes

PSF utilities in SimCADO
~~~~~~~~~~~~~~~~~~~~~~~~

In SimCADO we have pre-packaged some simulated PSFs which match the expectations of the 
different AO modes of MICADO. As new simulations of the AO capabilities become available we 
will include these new PSFs in SimCADO. The available PSFs are

- PSF_SCAO.fits
- PSF_MCAO.fits

which are the AO modes available. Additionally we provide a LTAO PSF (PSF_LTAO.fits) and a EELT diffraction limited PSF (PSF_POPPY.fits) calculated with `poppy <https://webbpsf.readthedocs.io/en/stable/poppy.html>`_ 

For the moment, the PSF is assumed constant accross the field. The new (refractored) version will
be capable to deal with field varying PSFs in a realistic manner. 

SimCADO also provides utility functions that are able to produce analytic PSFs to be used in the simulation in adition to the pre-calculated PSFs. Please check the :mod:`simcado.psf` module. The 
most important functions are the following:

- :func:`simcado.psf.poppy_eelt_psf` creates a diffraction limited PSF based with mirror segments provided by :func:`simcado.psf.get_eelt_segments`

- :func:`simcado.psf.seeing_psf` creates a seeing limited PSF with an user provided FWHM. Moffat and Gaussian profiles are available

- :func:`simcado.psf.poppy_ao_psf` creates an analytical AO PSF with a user provided Strehl ratio

All these functions can save the computed PSFs in fits format by specifying a filename. That file can be later used in the simulations as a parameter in :func:`simcado.simulation.run` using the ``SCOPE_PSF_FILE=filename`` keyword.

In the :doc:`notebook <6_Notebooks>` section you can find a few detailed examples how to create and work with these PSFs. 


Saving and reusing commands
---------------------------

The :class:`~.UserCommands` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Passing more than a few keyword-value pairs to the ``simcado.run()``
becomes tedious. SimCADO therefore provides a dictionary of commands so
that you can keep track of everthing that is happening in a simulation. ::

    >>> my_cmds = simcado.UserCommands()
    >>> simcado.run(my_src, cmds=my_cmds)

When initialised the :class:`~.UserCommands` object contains all the default
values for MICADO, as given in :doc:`Keywords <Z_Keywords>`. The
:class:`~.UserCommands` object is used just like a normal python dictionary:

::

    >>> my_cmds["OBS_EXPTIME"] = 180
    >>> my_cmds["OBS_EXPTIME"]
    180.0

It can be saved to disk and re-read later on:

::

    >>> my_cmds.writeto("path/to/new_cmds.txt")
    >>> new_cmds = simcado.UserCommands("path/to/new_cmds.txt")
    >>> new_cmds["OBS_EXPTIME"]
    180.0

If you prefer not to use interactive python and just want to dump a
commands file to edit in your favourite text editor:

::

    >>> simcado.commands.dump_defaults("path/to/cmds_file.txt")

More information on the :class:`~.UserCommands` object is given in the
:doc:`Examples Section <B_UserCommands>`

Behind the scenes of SimCADO
----------------------------

SimCADO uses 4 main classes during a simulation:

-  :class:`.Source` holds spatial and spectral information about the
   astronomical source of photons, e.g. galaxy, star cluster, etc.
-  ``OpticalTrain`` contains information on the various elements along
   the optical path, e.g. mirrors reflectivity curves, PSFs,
   instrumental distortion, etc.
-  ``Detector`` represents the focal plane detector array and contains
   information on the electronic characteristics of the detector chips
   and their physical positions.
-  :class:`~.UserCommands` is a dictionary of all the important keywords needed
   by SimCADO to run the simultationm, e.g. ``OBS_EXPTIME`` (exposure
   time) or ``INST_FILTER_TC`` (filter curve)

For more information on how SimCADO works please see the :doc:`SimCADO in Depth <A_DeepStuff>` section.

Things to watch out for
-----------------------

This space. It will soon expand!
