Examples with the Source object
===============================

The :class:`.Source` class is probably the most important class for testing
science cases. Therefore spending time on creating accurate :class:`.Source`
representations of the object of interest is key to getting good results
with SimCADO.

Basically a :class:`.Source` object represents photon sources using lists of
positions (``.x, .y``), a list of unique spectra (``.spectra``) and a
list of references which match each photon source to a spectrum in the
list of spectra (``.ref``). All sources (extended and point source) can
be decomposed into these lists. The advantage of using this approach is
that objects with highly similar spectra can both reference the same
position in ``.spectra``, thereby reducing the number of spectra that
need to be manipulated during a simulation.

.. contents::


My first Source object
----------------------

To begin with it is probably easiest to let SimCADO generate a
:class:`.Source` object. The convenience function ``simcado.source.star()``
will generate a :class:`.Source` object containing a single star. In this
case, we’ll choose an G2V star with a K-band magnitude of 20, placed 5
arcsec above the centre of the focal plane:

::

    >>> star_1 = simcado.source.star(mag=20, filter_name="K", spec_type="G2V", x=5, y=0))

``star_1`` the coordinates of the star are held in the arrays ``.x`` and
``.y``, and the G2V spectrum in \`.spectra´.

::

    >>> star_1.x[0], star_1.y[0]
    (5, 0)

The spectrum for ``star_1`` is held in ``.spectra`` and the central
wavelength of each of the spectral bins is in ``.lam``.

::

    >>> star_1.lam, star_1.spectra[0]
    <insert output here>

.. note::
    `.lam` is a (1,n) array where as `.spectra` is a (m,n) array where n is the number of bins in the spectra and m is the number of unique spectra in the `Source` object. If the `Source` only contains a single unique spectrum, then `.spectra` will be a (1,n) array too.

Combining :class:`.Source` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because a :class:`.Source` object is just a collection of arrays, it is easy
to add many together with the ``+`` operator:

::

    >>> star_2 = simcado.source.star(mag=22, filter_name="K", spec_type="A0V", x=2, y=-2)
    >>> two_stars = star_1 + star_2
    >>> two_stars.x, two_stars.y
    ((5, 2), (0, -2))

.. note::
    this is a very trivial example (for which SimCADO has a more elegant function: `simcado.source.stars()`), but it serves to illustrate the main idea. The overloaded `+` operator is very useful for combining objects, e.g. forground stars and background galaxies, in order to get a better representation of the sky.

.. note::
    if you are planning on creating sources for large numbers of stars (e.g. >>10 stars), using the plural function ``stars()`` will save you a lot of time.

Point sources
~~~~~~~~~~~~~

For generating a field of stars, SimCADO offers a series of convenience
functions. Please see the docstring or API documentation for more
information on how best to use them.

-  ``simcado.source.star()``
-  ``simcado.source.stars()``
-  ``simcado.source.star_grid()``
-  ``simcado.source.cluster()``

SimCADO’s in-built example spectra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Add a section on this

-  Stellar templates from Pickles (1998)

   -  ``simcado.source.SED()``
   -  ``simcado.source.empty_sky()``

-  Galaxy templates from STSCI

Each of these functions returns two arrays: ``lam`` and ``spec``

Using an image as a template for a :class:`.Source` object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If we have an extended source that we wish to simulate, e.g. a galaxy, a
nebula, etc. we can use the function
``simcado.source.source_from_image()``. The image must be a 2D
``numpy.ndarray``, but it can come from anywhere, e.g. a FITS file,
generated my another function, or even a MS Paint bitmap image.

::

    >>> image_1 = astropy.io.fits.getdata("orion.fits")
    >>> image_1
    array([[0.0, ... 0.0],
           ... ,
           [0.0, ... 0.0]])

Here SimCADO takes a the pixel coordinates of the image and converts
them to positions on the focal plane.

.. note::
    the user must specify a plate-scale in arcseconds (`pix_res=`) for the image. Each pixel with a value above a certain threshold (default `flux_threshold=0`) will be used in the `Source` object. The coordinates of these pixels are added to the arrays `.x` and `.y`. 

We also need to provide a spectrum for the image. This spectrum is
assumed to be the only spectrum for each pixel in the image. The pixel
values are then the intensity assigned to that spectrum at that pixel
position.

SimCADO provides the pickles library for stellar spectra. Unfortunately
there aren’t any built-in galactic spectra yet - for this the user will
need to provide their own spectrum.

::

    >>> lam, spec_1 = simcado.source.SED("G2V", "K", magnitude=20)
    >>> lam, spec_1
    (array([0.7 ... 2.5]), array([0.0 ... 0.0]))

With ``image_1``, ``lam`` and ``spec_1`` we can now build a :class:`.Source`
object for an orion-like nebula that has the spectrum of a sun-like
star.

::

    >>> simcado.source.source_from_image(image_1, lam, spec_1, pix_res=0.004, flux_threshold=0)

While this example is physically unrealistic, it serves the purpose of
showing how to build a :class:`.Source` object from an image. The user is

Images with multipe spectra
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In reality assigning a single spectrum to an extended object is of
limited use. For a :class:`.Source` to be realistic is should contain multiple
spectra for objects in different locations. The best way to simulate
this with SimCADO is to create a :class:`.Source` object for each unique group
of objects (e.g. old stellar population, star forming regions, AGN, etc)
and then combine them into a single :class:`.Source` object with the ``+``
operator.

As a worked example, lets create a “first-order” approximation to a star
forming galaxy. The two major components of this source are 1. the aged
stellar population and, 2. the star forming regions.

In our (very) crude model the aged stellar population can be approxiated
by an ellipse with Gaussian light distribution. As M stars make up the
majority of this population, we can assign a M0V spectrum to this
population.

::

    >>> from astropy.convolution import Gaussian2DKernel
    >>> from simcado.source import SED
    >>> 
    >>> old_pop = Gaussian2DKernel(128).array[::3,:]
    >>> m0v_spec = SED()

To illustrate (very crudely) the star forming regions we can create a
random distribution of elliptical Gaussians using the ``astropy``
function ``Gasussian2DKernel``:

Creating a :class:`.Source` object from scratch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To create a :class:`.Source` object from scratch, we initialise the object by
passing 5 (or 6) arrays. All the parameter names must be specified.

``sim.Source(lam=, spectra=, x=, y=, ref=, [weight=])``

where: + ``x, y`` - [each a ``numpy.ndarray``]. Coordinates for each
point source in the image in units of [arcsec] from the focal plane
centre

-  ``lam`` - [``numpy.ndarray``]. An array with the centre of the
   wavelength bins in [um] for each unique spectrum

-  ``spectra`` - [``numpy.ndarray``]. An (n, m) array holding n spectra,
   each with m values. Default units are [ph/s] Note - ``lam`` and
   ``spectra`` should use a constant bin width. Variable bin widths
   leads to unpredictable results.

-  ``ref`` - [``numpy.ndarray``]. An array to connect the point source
   at ``x[i]``, ``y[i]`` to a unique spectrum at ``spectra[j]``, i.e.
   ``ref[i] = j``

Optional keywords can be specified:

-  ``weight`` - [``numpy.ndarray``], optional. If two sources share the
   same spectrum, but are at different distances or have different
   luminosities a scaling factor can be specified to the spectrum when
   applied to each specific point source.
-  ``units`` [default ``"ph/s"``] is the units for the spectra, i.e. n
   phontons per second per spectral bin. The size of the spectral bins
   is resolution of the ``.lam`` array.

Combining two (or more) :class:`.Source` objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`.Source` objects can be created in different ways, but the underlying
table-structure is the same. Therefore adding :class:`.Source` objects
together means simply combining tables. The mathematical operator ``+``
can be used to do this:

::

    >>> # ... create a A0V star at (0,0) and a G2V star at (5,-5)
    >>> star_A0V = sim.source.star(20, spec_type="A0V", x=0, y=0)
    >>> star_G2V = sim.source.star(20, spec_type="G2V", x=5, y=-5)
    >>> 
    >>> src_combi = star_A0V + star_G2V
    >>> 
    >>> print(src_combi.x, src_combi.y)
    [0 5] [ 0 -5]

By adding different :class:`.Source` objects together, it is possible to build
up complex objects that will be representative of the observed sky,
e.g. old + new galaxy stellar population + gas emission + foreground
stars

See `examples <examples/Source>`__ for how to use the ``*`` and ``-``
operators with a :class:`.Source` object

Saving a :class:`.Source` object to disk
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`.Source` object is saved as a FITS file with two extensions. See
`How SimCADO works <in_depth/SimCADO>`__ for more on the file structure.

::

    >>> src_combi.write("my_src.fits")

The file can be read in at a later time by specifying ``filename=`` when
initialising a :class:`.Source` object - as stated above

::

    >>> my_src = sim.Source(filename="my_src.fits")
    

In-built :class:`.Source` object for a star cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a test object, SimCADO provides the function, with all distances in
parsecs
::

    sim.source.cluster(mass=1E4, distance=50000, half_light_radius=1)
    

SimCADO convenience functions
-----------------------------

* :func:`simcado.source.empty_sky`
* :func:`simcado.source.stars`
* :func:`simcado.source.cluster`
* :func:`simcado.source.SED`
* :func:`simcado.source.source_from_image`
* :class:`simcado.source.Source`
