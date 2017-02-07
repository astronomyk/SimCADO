Welcome to SimCADocs
====================

The (slowly expanding) documentation base for SimCADO

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
use SimCADO are available `here in the Notebooks section`_.

Reference Material
~~~~~~~~~~~~~~~~~~

-  The inner workings of SimCADO are described in detail in `Leschinski
   et al. (2016)`_

-  The current status of MICADO is described in `Davies et al. (2016)`_

Downloading and Installing
--------------------------

For more information, see the `Downloads`_ section

**SimCADO has only been tested in Python 3.x**.

It is hightly recommended to use Python 3, however the basics of
generating images will still work in Python 2.7. We cannot guarantee
this though. See the `Features`_ page for more info on which functions
with which Python version.

The quick way:

::

    $ pip3 install --user http://www.univie.ac.at/simcado/SimCADO.tar.gz

The **first time** in python

::

    >>> import simcado
    >>> simcado.get_extras()
    >>>
    >>> # !! Only works in Python 3 - See Downloads section
    >>> simcado.install_noise_cube()

Keeping SimCADO updated
~~~~~~~~~~~~~~~~~~~~~~~

As MICADO developes, the data files that SimCADO uses will also be
updated. Therefore before you do any major work with SimCADO we *HIGHLY*
recommend calling:

::

    >>> simcado.get_extras()

Running a simulation in 3 lines
-------------------------------

The easiest way to run a simulation is to create, or load, a Source
object and then call the ``.run()`` command. If you specify a filename,
the resulting image will be output to a FITS file under that name. If
you do not specify a filename, the output will be returned to the
console/notebook as an ``astropy.io.fits.HDUList`` object.

To begin, we will import the simcado module (assuming it is already
installed).

::

    >>> import simcado

At the very least, we need to create a ``Source`` object which contains
both spatial and spectral information on our object of interest. Here we
use the built-in command ``.source.cluster()`` to create a
``Source``-object for a 10000-Msun stellar cluster. (See `Creating
Sources`_ for more information).

::

    >>> src = simcado.source.cluster()

We now pass the ``source`` object through SimCADO. This is as easy as
calling ``.run()``. If we specify a ``filename``, SimCADO

.. _here in the Notebooks section: examples/Notebooks.md
.. _Leschinski et al. (2016): https://arxiv.org/pdf/1609.01480v1.pdf
.. _Davies et al. (2016): https://arxiv.org/pdf/1607.01954.pdf
.. _Downloads: Download.md
.. _Features: Features.md
.. _Creating Sources: examples/Source.md