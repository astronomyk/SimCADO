.. SimCADO documentation master file, created by
   sphinx-quickstart on Wed Nov 21 15:05:49 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SimCADocs!
=====================

The (slowly expanding) documentation base for SimCADO, the instrument data
simulation package for MICADO at the ELT.

.. figure:: _static/images/Pulsar_light_curve_Davies.png
    :figwidth: 600 px
    :align: center

    This is the lightcurve that follows the shape of the Crab Pulsar, but
    scaled so that each pulse is easily visible in a 5ms windowed H-band
    exposure with MICADO. The star is H=15mag. The light curve has been
    extracted from the intermediate non-destructive reads from an up-the-ramp
    simulated exposure.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   user_docs/1_GettingStarted
   user_docs/2_Download
   user_docs/3_FAQs
   user_docs/4_Updates
   user_docs/5_Gallery
   user_docs/6_Notebooks

   user_docs/A_DeepStuff
   user_docs/B_UserCommands
   user_docs/C_Source
   user_docs/Z_Keywords

   reference/simcado

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. All Notebooks and API available on the University of Vienna's server:
   www.univie.ac.at/simcado

SimCADO in a nutshell
---------------------

SimCADO is a python package designed to simulate the effects of the
Atmosphere, E-ELT, and MICADO instrument on incoming light. The current
version (v0.5) can simulate the MICADO imaging modi (4mas and 1.5mas per
pixel in the wavelength range 0.7µm to 2.5µm).

iPython/Jupyter notebooks
~~~~~~~~~~~~~~~~~~~~~~~~~

A (continualy expanding) series of iPython Notebooks detailing how to
use SimCADO are available :doc:`here in the Notebooks
<user_docs/6_Notebooks>`
section.

.. hint::
    Don't feel like sifting through documentation? Common commands and examples
    are on the SimCADO cheat-sheet:

    * `PDF version`_ or
    * `Jupyter Notebook`_
    * `Presentation for October 2017 Science Team Meeting`_


.. _PDF version: ./_static/pdfs/SimCADO_cheatsheet.pdf
.. _Jupyter Notebook: https://nbviewer.jupyter.org/github/astronomyk/SimCADO/blob/master/docs/source/_static/python_notebooks/SimCADO-cheat-sheet.ipynb
.. _Presentation for October 2017 Science Team Meeting: ./_static/pdfs/SimCADO_status_Oct_2017.pdf

Reference Material
~~~~~~~~~~~~~~~~~~

-  The inner workings of SimCADO are described in detail in `Leschinski
   et al. (2016)`_

.. _Leschinski et al. (2016): https://arxiv.org/pdf/1609.01480v1.pdf

-  The current status of MICADO is described in `Davies et al. (2016)`_

.. _Davies et al. (2016): https://arxiv.org/pdf/1607.01954.pdf

Downloading and Installing
--------------------------

For more information, see the :doc:`Download <user_docs/2_Download>` section

**SimCADO has only been tested in Python 3.x**.

It is hightly recommended to use Python 3, however the basics of
generating images will still work in Python 2.7. We cannot guarantee
this though.

The quick way::

    $ pip3 install --user http://www.univie.ac.at/simcado/SimCADO.tar.gz

The **first time** in python you should update the SimCADO data directory with
:func:`~.simcado.utils.get_extras`::

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
recommend calling :func:`~.simcado.utils.get_extras`::

    >>> simcado.get_extras()

Running a simulation in 3 lines
-------------------------------

The easiest way to run a simulation is to create, or load, a Source
object and then call the :func:`.run` command. If you specify a filename,
the resulting image will be output to a FITS file under that name. If
you do not specify a filename, the output will be returned to the
console/notebook as an :class:`~.astropy.io.fits.hdu.hdulist.HDUList` object.

To begin, we will import the simcado module (assuming it is already
installed).::

    >>> import simcado

At the very least, we need to create a :class:`.Source` object which contains
both spatial and spectral information on our object of interest. Here we
use the built-in command :func:`simcado.source.cluster()` to create a
:class:`.Source` object for a 10000-Msun stellar cluster. (:doc:`Creating
Sources <user_docs/C_Source>` for more information).::

    >>> src = simcado.source.cluster()

We now pass the :class:`.Source` object through SimCADO. This is as easy as
calling :func:`.run`. If we specify a ``filename``, SimCADO will write the
output to disk in the form of a FITS file. If no ``filename`` is given, then
SimCADO returns an astropy :mod:`~.astropy.io.fits` object to the console or
notebook.::

    >>> simcado.run(src, filename="my_first_sim.fits")


Changing simulation parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :func:`.run` also takes any :doc:`configuration keywords <user_docs/Z_Keywords>` as
parameters for running the simulation. For example, the default exposure time
for the simulation is 60 seconds, however this can be increased of decreased by
using the keyword `OBS_EXPTIME` (and/or combining it with `OBS_NDIT`). A stacked
6x 10 minute observation sequence would look like::

    >>> simcado.run(src, filename="my_first_sim.fits", OBS_EXPTIME=600, OBS_NDIT=6)

That's it. Of course SimCADO can also go in the other direction, providing many
more levels of complexity, but for that the reader is directed to the examples
pages and/or the :doc:`API documentation <reference/simcado>`.

SimCADO building blocks
-----------------------
For a brief explanation of how SimCADO works and which classes are relevant,
please see either the :doc:`GettingStarted <user_docs/1_GettingStarted>` or
:doc:`SimCADO in depth <user_docs/A_DeepStuff>` section.

Bugs and Issues
---------------

We freely admit that there may still be several bugs that we have not found.
If you come across an buggy part of SimCADO, *please please* tell us. We can't
make SimCADO better if we don't know about things.

The preferable option is to open an issue on our Github page:
`gastronomyk/SimCADO/issues`_, or you can contact either one of us directly.

.. _gastronomyk/SimCADO/issues: https://github.com/gastronomyk/SimCADO/issues,

Contact
-------

For questions and complaints alike, please contact the authors:

* kieran.leschinski@univie.ac.at
* oliver.czoske@univie.ac.at

**Developers (Vienna):** Kieran Leschinski, Oliver Czoske

**Data Flow Team Leader (Gronigen):** Gijs Verdoes Kleijn

**MICADO home office (MPE):** http://www.mpe.mpg.de/ir/micado
