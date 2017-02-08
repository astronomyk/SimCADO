Welcome to SimCADocs
====================

The (slowly expanding) documentation base for SimCADO
:mod:`astropy.io`

.. figure:: images/Omega_Cen_Fabricius.png
    :figwidth: 600 px
    :align: center

    Omega Cen as imaged with HST/WFC3, HST/SimCADO and MICADO/SimCADO by
    Maximilian Fabricius (MPE). The synthetic images
    of the same region of Omega Cen are based on the HST catalog by
    Anderson & van der Marel 2010 and augmented by all the faint stars
    that did not end up in the HST catalogue.


SimCADO in a nutshell
---------------------

SimCADO is a python package designed to simulate the effects of the
Atmosphere, E-ELT, and MICADO instrument on incoming light. The current
version (v0.2) can simulate the MICADO imaging modi (4mas and 1.5mas per
pixel in the wavelength range 0.7µm to 2.5µm).

iPython/Jupyter notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~

A (continualy expanding) series of iPython Notebooks detailing how to
use SimCADO are available :doc:`here in the Notebooks <examples/Notebooks>` section.

Reference Material
~~~~~~~~~~~~~~~~~~

-  The inner workings of SimCADO are described in detail in `Leschinski
   et al. (2016)`_

-  The current status of MICADO is described in `Davies et al. (2016)`_

Downloading and Installing
--------------------------

For more information, see the :doc:Downloads section

**SimCADO has only been tested in Python 3.x**.

It is hightly recommended to use Python 3, however the basics of
generating images will still work in Python 2.7. We cannot guarantee
this though. See the `Features`_ page for more info on which functions
with which Python version.

The quick way:

::

    $ pip3 install --user http://www.univie.ac.at/simcado/SimCADO.tar.gz

The **first time** in python you should update the SimCADO data directory with
:func:`~.simcado.utils.get_extras`:

::

    >>> import simcado
    >>> simcado.get_extras()
    >>>
    >>> # !! Only works in Python 3 - See Downloads section
    >>> simcado.install_noise_cube()

If you running Python 3, it would be helpful to expand the internal detector 
noise cube with :func:`.install_noise_cube`. 

    
Keeping SimCADO updated
~~~~~~~~~~~~~~~~~~~~~~~

As MICADO developes, the data files that SimCADO uses will also be
updated. Therefore before you do any major work with SimCADO we **HIGHLY**
recommend calling :func:`~.simcado.utils.get_extras`:

::

    >>> simcado.get_extras()

Running a simulation in 3 lines
-------------------------------

The easiest way to run a simulation is to create, or load, a Source
object and then call the :func:`.run` command. If you specify a filename,
the resulting image will be output to a FITS file under that name. If
you do not specify a filename, the output will be returned to the
console/notebook as an :class:`~.astropy.io.fits.hdu.hdulist.HDUList` object.

To begin, we will import the simcado module (assuming it is already
installed).
::

    >>> import simcado

At the very least, we need to create a :class:`.Source` object which contains
both spatial and spectral information on our object of interest. Here we
use the built-in command :func:`simcado.source.cluster()` to create a
:class:`.Source`-object for a 10000-Msun stellar cluster. (See `Creating
Sources`_ for more information).
::

    >>> src = simcado.source.cluster()

We now pass the :class:`.Source` object through SimCADO. This is as easy as
calling :func:`.run`. If we specify a ``filename``, SimCADO will write the 
output to disk in the form of a FITS file. If no ``filename`` is given, then 
SimCADO returns an astropy :mod:`~.astropy.io.fits` object to the console or 
notebook.
::

    >>> simcado.run(src, filename="my_first_sim.fits")

    
Changing simulation parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :func:`.run` also takes any `configuration keywords`_ as parameters for running the simulation. For example, the default exposure time for the simulation is 60 seconds, however this can be increased of decreased by using the keyword `OBS_EXPTIME` (and/or combining it with `OBS_NDIT`). A stacked 6x 10 minute observation sequence would look like:

    >>> simcado.run(src, filename="my_first_sim.fits", OBS_EXPTIME=600, OBS_NDIT=6)
    
That's it. Of course SimCADO can also go in the other direction, providing many more levels of complexity, but for that the reader is directed to the examples pages and/or the [API](API/_build/index.html) documentation

SimCADO building blocks
-----------------------
For a brief explanation of how SimCADO works and which classes are relevant, please see either the [Getting Started](GettingStarted.md) or [SimCADO in depth](deep_stuff/SimCADO.md) section.

Bugs and Issues
---------------

We freely admit that there may still be several bugs that we have not found. If you come across an buggy part of SimCADO, *please please* tell us. We can't make SimCADO better if we don't know about things.

The preferable option is to open an issue on our Github page: `gastronomyk/SimCADO/issues`_, or you can contact either one of us directly.

Contact
-------

For questions and complaints alike, please contact the authors:

* kieran.leschinski@univie.ac.at
* oliver.czoske@univie.ac.at

**Developers (Vienna):** Kieran Leschinski, Oliver Czoske

**Data Flow Team Leader (Gronigen):** Gijs Verdoes Kleijn

**MICADO home office (MPE):** http://www.mpe.mpg.de/ir/micado




.. _Leschinski et al. (2016): https://arxiv.org/pdf/1609.01480v1.pdf
.. _Davies et al. (2016): https://arxiv.org/pdf/1607.01954.pdf
.. _Downloads: Download.md
.. _Features: Features.md
.. _Creating Sources: examples/Source.md
.. _configuration keywords: Keywords.md
.. _gastronomyk/SimCADO/issues: https://github.com/gastronomyk/SimCADO/issues,